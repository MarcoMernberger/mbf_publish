from pathlib import Path
import pypipegraph as ppg
from mbf_externals.util import write_md5_sum
import pickle
import json


def prep_scb(*objects_to_push):
    """Main entry point to publish data to our scb system
    """
    objs = []
    for x in objects_to_push:
        if isinstance(x, dict):
            objs.extend(x.values())
        elif hasattr(x, '__iter__'):
            objs.extend(iter(x))
        else:
            print(x)
            print(dir(x))
            raise ValueError()
    submission = SCBSubmission(objs)
    return submission.dump_meta_data()


class SCBSubmission:
    def __init__(self, objs):
        self.output_path = Path('web/scb')
        self.output_path.mkdir(exist_ok=True, parents=True)
        self.objs = objs
        self.used = set()
        self.genomes = set()
        self.meta_data = {}
        self.deps = []
        self.errors = []
        self.extract_genomes()
        self.extract_lanes()
        if len(self.used) != len(self.objs):
            missing = set(self.objs) - self.used
            raise ValueError("unused objects", missing)
        if self.errors:
            raise ValueError(self.errors)

    def extract_genomes(self):
        for o in self.objs:
            if hasattr(o, "genome"):
                self.genomes.add(o)

    def dump_meta_data(self):
        self.dump_meta_data_json()
        self.dump_rsync_list()

    def dump_meta_data_json(self):
        output_filename = "web/scb/metadata.json"

        def descend_and_replace_callbacks(d):
            print('descend_and_replace_callbacks')
            for k,v in d.items():
                if hasattr(v, '__call__'):
                    print('bingo')
                    d[k] = v()
                elif isinstance(v, dict):
                    descend_and_replace_callbacks(v)
                elif isinstance(v, list):
                    for vx in v:
                        if isinstance(vx, dict):
                            descend_and_replace_callbacks(vx)

        def dump(output_filename):
            descend_and_replace_callbacks(self.meta_data)
            with open(output_filename, "w") as op:
                json.dump(self.meta_data, op, indent=4)  # make sure it's sane, no funny objects

        return ppg.FileGeneratingJob(output_filename, dump).depends_on(self.deps)

    def dump_rsync_list(self):
        output_filename = 'web/scb/rsync_list.txt'
        def dump(output_filename):
            paths_to_copy = set()
            for group in self.meta_data:
                for entry in self.meta_data[group]:
                    for x in ['table_path', 'path_bigbed', 'path_table', 'path_bam', 'path']:
                        if x in entry:  # genes
                            paths_to_copy.add(Path(entry[x]).absolute())
                    if 'path_bam' in entry:
                            paths_to_copy.add(Path(entry['path_bam'] + '.bai').absolute())
            for fn in Path('web/scb').glob("*"):
                paths_to_copy.add(fn.absolute())
            output = ''
            paths_to_copy.add("/project/web/scb/metadata.json")
            for x in sorted([str(x) for x in paths_to_copy]):
                if x.startswith('/project'):
                    output += x[len('/project/'):]
                    output += "\n"
            Path(output_filename).write_text(output)


        return ppg.FileGeneratingJob(output_filename, dump).depends_on(self.deps)

    def extract_lanes(self):
        self.meta_data["lanes"] = []
        for entry in self.objs:
            if hasattr(entry, "get_bam_names"):
                vids = flatten_vid(entry.vid, entry, self.errors)
                fn = entry.get_bam_names()[0]
                md5 = get_md5(fn, entry.load())
                self.meta_data["lanes"].append(
                    {
                        "vids": vids,
                        "path_bam": fn,
                        "displayname": get_scb_name(entry),
                        "md5": lambda: Path(md5[1]).read_text(),
                        "md5_path": fn,
                        "scb_comment": (
                            entry.scb_comment if hasattr(entry, "scb_comment") else ""
                        ),
                        "id": (
                            entry.scb_id if hasattr(entry, "scb_id") else None
                        ),  # this allows you to permanently tie it to an entry in SCBs databases
                        #"browser_options": entry.gbrowse_options,
                        #"has_splicing": entry.has_spliced_reads,
                        "genome_label": extract_genome_label(entry),
                    }
                )
                self.deps.append(md5[0])
                self.register_used(entry)



    def register_used(self, entry):
        if entry in self.used:
            raise ValueError("double use", entry)
        self.used.add(entry)


def get_md5(input_name, input_job):
    input_name = Path(input_name)
    target_name = input_name.with_name(input_name.name + ".md5sum")

    def gen(target_filename):
        target_filename = Path(target_filename)
        write_md5_sum(target_filename.with_name(target_filename.stem))

    return (
        ppg.FileGeneratingJob(target_name, gen).depends_on(input_job),
        str(target_name),
    )


def flatten_vid(vid, entry, errors, collector=None):
    if collector is None:
        collector = []
    if vid is None:
        errors.append(
            (
                entry.name if hasattr(entry, "name") else entry.column_name,
                "Vid was none - all samples have to be associated with a vid",
            )
        )
    elif hasattr(vid, "__iter__") and not isinstance(vid, str):
        for v in vid:
            try:
                flatten_vid(v, entry, errors, collector)
            except RecursionError:
                print(entry)
                raise
    else:
        collector.append(vid)
    return collector


def get_scb_name(obj):
    if hasattr(obj, "scb_name") and obj.scb_name:
        return obj.scb_name
    elif hasattr(obj, "sheet_name") and obj.sheet_name:
        return obj.sheet_name + "/" + obj.name
    else:
        return obj.name


def extract_genome_label(obj):
    if obj.genome.species == "Mus_musculus":
        return "mm10"
    elif obj.genome.species == "Homo_sapiens":
        if int(obj.genome.revision) >= 77:
            return "hg38"
        else:
            return "hg19"
    elif obj.genome.species == "Drosophila_melanogaster":
        return "dm89"
    elif obj.genome.species == "Ustilago_maydis":
        return "um33"
    else:
        raise ValueError(
            "don't know how to tell the scb about this genome: %s" % obj.genome
        )
