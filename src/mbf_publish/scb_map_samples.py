import requests
import sys
import os
import itertools
import re
from pathlib import Path
import json
from collections import OrderedDict, namedtuple
import tabulate

SampleInfo = namedtuple("SampleInfo", ["name", "fastqs", "vid"])


def map_samples(
    sample_paths,
    scb_project_ids,
    mapping_filename=".scb_sample_matches",
    filter_fastqs=None,
    expect_vids_missing=None,
):
    """(interactively) map fastq.gz samples and vids from the scb.
    Once it's done, return a dictionary:
    {sample_name -> {
        'vid': 'XYZ000',
        'fastqs': [path, path, path...],
        'display_name': 'sample_name'
    }
    The result is stored in @mapping_file and reused the next time
    (without interactivity if complete).
    You can force an interactive run by opening @mapping_file and changing
    force_rebuild to be true (json format).

    If you expect vids to be missing,
    you can add as a set in @expect_vids_missing,
    and they're filtered from the scb response

    If you have discoverable fastqs that you want to ignore,
    pass their 'names' (ie. the part before _S\\d+_L\\d+)
    in @filter_fastqs.
    """
    sample_paths = [Path(x) for x in sample_paths]
    mapping_filename = Path(mapping_filename)

    fastq_samples = discover_fastq_samples(sample_paths)
    if filter_fastqs:
        for k in filter_fastqs:
            if k in fastq_samples:
                print(f"Not using fastq: {k}")
                del fastq_samples[k]
    if not fastq_samples:
        raise ValueError("No fastqs found in {sample_paths}")

    existing_mapping, require_rebuild = load_mapping(mapping_filename, fastq_samples)
    if "--map-samples" in sys.argv:
        require_rebuild = True
        default_hide = True
    else:
        default_hide = False

    if not require_rebuild and existing_mapping_complete(
        existing_mapping, fastq_samples
    ):
        result = existing_mapping
    else:
        vid_information = fetch_vid_info(scb_project_ids)
        if expect_vids_missing:
            for vid in expect_vids_missing:
                if vid in vid_information:
                    print("Filtered", vid, vid_information[vid])
                    del vid_information[vid]
        result = map_interactive(
            existing_mapping,
            fastq_samples,
            vid_information,
            mapping_filename,
            default_hide,
        )
    return [
        SampleInfo(name=v["display_name"], vid=v["vid"], fastqs=fastq_samples[k])
        for (k, v) in result.items()
    ]


def existing_mapping_complete(existing_mapping, fastq_samples, print_errors=False):
    if existing_mapping.keys() != fastq_samples.keys():
        if print_errors:
            missing = set(fastq_samples.keys()).difference(existing_mapping.keys())
            print(f"Not all fastqs associated. Missing: {missing}")
        return False
    else:
        for k, value in existing_mapping.items():
            # these_fastqs = value.get("fastqs", False)
            # if these_fastqs is False:
            # print("Mapping for non-existant fastq?!")
            # return False
            # if set(fastq_samples[k]) != set(these_fastqs):
            # return False
            if not "vid" in value:
                print(f"No vid assigned for {k}")
                return False
            if value.get("display_name", "").strip() == "":
                print(f"No display name set for {k}")
                return False
            if "/" in value["display_name"]:
                print(f"Invalid character in name for {k} - {value['display_name']}")
                return False
    return True


def discover_fastq_samples(sample_paths):
    r"""Discover fastq.gz files. Assume that anything before _S\d+_L\d+_ is the sample name"""
    by_key = {}
    for p in sample_paths:
        for candidate in p.glob("**/*.fastq.gz"):
            rx = r"_S\d+_L\d+_"
            match = re.search(rx, str(candidate.name))
            if not match:
                raise ValueError(f"Fastq -> sample name failed for {candidate}")
            name = candidate.name[: match.start()]
            if not name in by_key:
                by_key[name] = set()
            by_key[name].add(candidate)

    res = OrderedDict()
    for k in sorted(by_key):
        res[k] = by_key[k]
    return res


def load_mapping(path, fastq_samples):
    mapping = {}
    force_rebuild = True
    if path.exists():
        try:
            with open(path) as op:
                j = json.load(op)
            force_rebuild = bool(j.get("force_rebuild", True))
            mapping = j.get("mapping", {})
        except json.JSONDecodeError:
            print("could not read existing mapping file")
            pass
    for k in list(mapping.keys()):
        if not k in fastq_samples:
            del mapping[k]
    return mapping, force_rebuild


def save_mapping(mapping, path):
    with open(path, "w") as op:
        json.dump(
            {"mapping": mapping, "force_rebuild": False}, op, indent=4, sort_keys=True
        )


def fetch_vid_info(scb_project_ids):
    url = "https://mbf.imt.uni-marburg.de/scb/show_experiment/%i?machine_readable=True"
    auth = requests.auth.HTTPBasicAuth(
        os.environ["MBF_AUTH_USER"], os.environ["MBF_AUTH_PASSWORD"]
    )

    result = {}
    for id in scb_project_ids:
        r = requests.get(url % id, auth=auth)
        if r.status_code == 404:
            raise KeyError("Invalid experiment id: %i not found on server" % id)
        r = json.loads(r.text)
        for k, v in r["samples"].items():
            v["project"] = r["short_name"]
            v["project_id"] = id
            v["vid"] = k
            result[k] = v
    return result


def map_interactive(  # noqa:C901
    existing_mapping, fastq_samples, vid_information, mapping_filename, default_hide
):
    show_fastqs = not default_hide
    show_vids = not default_hide

    while True:
        mapped_vids = {v["vid"]: k for (k, v) in existing_mapping.items()}
        if existing_mapping_complete(existing_mapping, fastq_samples):
            print("Remapping request")
        else:
            print("Mapping incomplete")
        print("")

        print_mapping(existing_mapping)
        print("Vid infomation")
        if show_vids:
            print_vid_information(vid_information, mapped_vids)
        else:
            print("\tHidden - use togglevid/show to show")
        print("Fastqs")
        if show_fastqs:
            print_fastqs(fastq_samples, existing_mapping)
        else:
            print("\tHidden - use togglefastq/show to show")
        print("h<enter> for help")
        cmd = sys.stdin.readline().strip()
        map_match = re.match("([0-9]+)([A-Za-z]{2}[0-9]{1,3})(d|f|.+)?", cmd)
        if re.match("^-[A-Za-z]{2}[0-9]{1,3}$", cmd):
            # unmap
            vid = parse_vid(cmd[1:])
            for k, v in existing_mapping.items():
                if v.get("vid", "") == vid:
                    v["vid"] = ""
                    print(f"unmatched {vid} from {k}")
                    break
            else:
                print("could not find a mapping for {vid}")
        elif re.match("^-[0-9]", cmd):
            # unmap based on fastq
            fastq_no = int(cmd[1:])
            if fastq_no >= len(fastq_samples):
                print(f"invalid fastq number")
                continue
            key = list(fastq_samples)[fastq_no]
            if key in existing_mapping:
                existing_mapping[key]["vid"] = ""
        elif map_match:
            # map
            fastq_no = int(map_match.groups()[0])
            vid = parse_vid(map_match.groups()[1])
            if not vid in vid_information:
                print(f"invalid vid {vid}")
                continue
            if fastq_no >= len(fastq_samples):
                print(f"invalid fastq number")
                continue
            key = list(fastq_samples)[fastq_no]

            assign = False
            description = map_match.groups()[2]
            if description is None:
                description = "d"
            description = description.strip()
            if description == "d" or description == "":
                description = vid_information[vid]["user_description"]
                assign = True
            elif description == "f":
                description = key
                assign = True
            if not key in existing_mapping:
                existing_mapping[key] = {
                    "vid": vid,
                    "display_name": description,
                    "factors": {
                        k: vid_information[vid][k]
                        for k in extract_factors(vid_information[vid])
                    },
                }
            else:
                existing_mapping[key]["vid"] = vid
                existing_mapping[key]["factors"] = {
                    k: description for k in extract_factors(vid_information[vid])
                }
                if assign:
                    existing_mapping[key]["display_name"] = description
        elif re.match("[A-Za-z]{2}[0-9]{1,3}f", cmd):
            # name from vid and fastq name
            vid = parse_vid(cmd[:-1])
            change_display_name(existing_mapping, vid, lambda key, mapping: key)
        elif re.match("[A-Za-z]{2}[0-9]{1,3}d", cmd):
            # name from vid and user description
            vid = parse_vid(cmd[:-1])
            change_display_name(
                existing_mapping,
                vid,
                lambda key, mapping: vid_information[mapping["vid"]][
                    "user_description"
                ],
            )
        elif re.match("[A-Za-z]{2}[0-9]{1,3}.+", cmd):
            # name from vid and custom
            vid, text = re.match("([A-Za-z]{2}[0-9]{1,3})(.+)", cmd).groups()
            vid = parse_vid(vid)
            text = text.strip()
            change_display_name(existing_mapping, vid, lambda key, mapping: text)
        elif re.match("[0-9]+f", cmd):
            # name from vid and fastq name
            fastq_no = int(cmd[:-1])
            if fastq_no >= len(fastq_samples):
                print(f"invalid fastq number")
                continue
            key = list(fastq_samples)[fastq_no]
            try:
                vid = existing_mapping[key]["vid"]
                change_display_name(existing_mapping, vid, lambda key, mapping: key)
            except KeyError:
                print(f"no mapping for {fastq_no}")

        elif re.match("[0-9]+d", cmd):
            # name from vid and user description
            fastq_no = int(cmd[:-1])
            if fastq_no >= len(fastq_samples):
                print(f"invalid fastq number")
                continue
            key = list(fastq_samples)[fastq_no]
            try:
                vid = existing_mapping[key]["vid"]
                change_display_name(
                    existing_mapping,
                    vid,
                    lambda key, mapping: vid_information[mapping["vid"]][
                        "user_description"
                    ],
                )
            except KeyError:
                print(f"no mapping for {fastq_no}")

        elif re.match("[0-9]+.+", cmd):
            # name from vid and user description
            fastq_no, text = re.match("([0-9]+)(.+)", cmd).groups()
            fastq_no = int(fastq_no)
            text = text.strip()
            if fastq_no >= len(fastq_samples):
                print(f"invalid fastq number")
                continue
            key = list(fastq_samples)[fastq_no]
            try:
                vid = existing_mapping[key]["vid"]
                change_display_name(existing_mapping, vid, lambda key, mapping: text)
            except KeyError:
                print(f"no mapping for {fastq_no}")
        elif cmd.startswith("s/"):
            parts = cmd.replace("\\/", "%%SAVED_SLASH%%").split("/")
            parts = [x.replace("%%SAVED_SLASH%%", "/") for x in parts]
            if len(parts) != 2 and len(parts) != 3:
                print(
                    "invalid s/ command. Must be s/search/replace, escape internal / with \\"
                )
            else:
                if len(parts) == 2:
                    replace = ""
                else:
                    replace = parts[2]
                try:
                    search = parts[1]
                    search = search.replace("\\/", "/")
                    print("replacing %s with %s" % (repr(search), repr(replace)))
                    for v in existing_mapping.values():
                        v["display_name"] = re.sub(search, replace, v["display_name"])
                except re.error as e:
                    print(e)
                    continue
        elif cmd == "abort" or cmd == "q":
            print("aborting")
            sys.exit(0)
        elif cmd == "reset":
            print("resetting mapping")
            existing_mapping = {}
        elif cmd == "h" or cmd == "help":
            print_help()
        elif cmd == "togglevid":
            show_vids = not show_vids
        elif cmd == "togglefastq":
            show_fastqs = not show_fastqs
        elif cmd == "show":
            show_vids = True
            show_fastqs = True
        elif cmd == "hide":
            show_vids = False
            show_fastqs = False
        elif cmd == "done":
            if not existing_mapping_complete(existing_mapping, fastq_samples, True):
                print("Mapping not complete! Not leaving")
            else:
                break
        elif cmd == "":
            continue
        else:
            print("I did not understand this")
            continue
        save_mapping(existing_mapping, mapping_filename)
    return existing_mapping


def parse_vid(proto_vid):
    """Turn rm0 into RM000 etc"""
    letters = proto_vid[:2].upper()
    numbers = proto_vid[2:]
    numbers = "0" * (3 - len(numbers)) + numbers
    return letters + numbers


def change_display_name(existing_mapping, vid, name_callback):
    for k, v in existing_mapping.items():
        if v.get("vid", "") == vid:
            d = name_callback(k, v)
            if d:
                v["display_name"] = d
                print(f"assigned name {v['display_name'] } to {vid}")
            else:
                print("Not assigning an empty display name")
            break
        else:
            print("could not find a mapping for {vid}")


def print_mapping(mapping):
    print("Mapping")
    if not mapping:
        print("\t (no mapping)")
    else:
        header = ["Fastq", "Vid", "Name", "info"]
        entries = []
        for k in sorted(mapping, key=lambda k: mapping[k]["vid"]):
            entries.append(
                (
                    k,
                    mapping[k]["vid"],
                    mapping[k]["display_name"],
                    format_factors(mapping[k]["factors"]),
                )
            )
        print(tabulate.tabulate(entries, headers=header))
    print("")


def print_vid_information(vid_information, mapped_vids):
    for project_name, vids in itertools.groupby(
        vid_information.values(), lambda x: x["project"]
    ):
        vids = sorted(vids, key=lambda x: x["vid"])
        print(f"{project_name} - {vids[0]['project_id']}")
        # print("\t\tMapped\tVid\tInfo")
        entries = []
        header = ["fastq", "vid"]
        for vinfo in vids:
            if vinfo["vid"] in mapped_vids:
                v = mapped_vids[vinfo["vid"]]
            else:
                v = ""
            e = [v, vinfo["vid"]]
            for k in extract_factors(vinfo):
                e.append(vinfo[k])
            if len(header) == 2:
                for k in extract_factors(vinfo):
                    header.append(k)
            entries.append(e)
        print(indent(tabulate.tabulate(entries, headers=header), 1))
        print("")


def indent(text, count):
    lines = text.split("\n")
    return "\n".join(("\t" * count + l for l in lines))


def extract_factors(vid_info):
    keys = sorted(vid_info.keys())
    for k in ("vid", "project_id", "project"):
        if k in keys:
            keys.remove(k)
    return keys


def format_factors(vid_info):
    keys = extract_factors(vid_info)
    res = ""
    for k in keys:
        res += f"{k}={vid_info[k]}; "
    return res


from itertools import zip_longest


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def print_fastqs(fastqs, existing_mapping):
    import math

    entries = []
    for ii, key in enumerate(fastqs):
        vid = existing_mapping.get(key, {}).get("vid", "")
        key = str(ii), key, vid, ""
        entries.append(key)
    entries = list(grouper(entries, math.ceil(len(entries) / 3.0), ("", "", "", "")))
    entries = list(zip(*entries))
    entries = [list(itertools.chain(*e)) for e in entries]
    print(tabulate.tabulate(entries, headers=["no", "fastq", "vid", "  "] * 3))


def print_help():
    print("Basic command () = optional, confirm with enter")
    print('\t"(nn)XYmmm(d|f|description)"')
    print("\t\t Map fastq nn to vid XYmmm, ")
    print(
        "\t\toptionally assign user_description (default) or fastq name or full description"
    )
    print("\tUnmap: -XYmmm or -nn")
    print("\tSet display name:")
    print("\t\tto fastq: XYmmmf or nnf or nnd (nn must be assigned)")
    print("\t\tto user_description: XYmmmd or nnd (nn must be assigned)")
    print("\t\tarbitrary XYmmm <your text here> or nn <text here>")
    print("\tfull fastq info: f")
    print("\tabort (but keep mapping): abort or q")
    print("\treset mapping: reset")
    print("\toggle vid information visibility: togglevid")
    print("\toggle fastq visibility: togglefastq")
    print("\thide all: hide")
    print("\tshow all: show")
    print("\ts/search/replace (regexps search)")
    print("\tFinish: done - only leaves if all mappings are completed!")
    print("")
