from pathlib import Path
import pypipegraph as ppg
from mbf_externals.utils import write_md5_sum
import pickle
import json


def prep_scb(*objects_to_push):
    """Main entry point to publish data to our scb system
    """
    objs = []
    for x in objects_to_push:
        if isinstance(x, dict):
            objs.extend(x.values())
        if hasattr(x, iter):
            objs.extend(iter(x))
    submission = SCBSubmission(objs)
    return submission.dump_metadata()


class SCBSubmission:
    def __init__(self, objs):
        self.output_path = Path('web/scb')
        self.output_path.mkdir(exist_ok=True, parents=True)
        self.objs = objs
        self.genomes = set()
        self.extract_genomes()
        self.meta_data = {}
        self.deps = []
        self.extract_lanes(objs)
        self.dump_meta_data()

    def extract_genomes(self):
        for o in self.objs:
            if hasattr(o, "genome"):
                self.genomes.add(o)

    def dump_meta_data(self):
        output_filename = "web/scb/metaddata.dat"

        def dump(output_filename):
            json.dumps(self.meta_data)  # make sure it's sane, no funny objects
            with open(output_filename, "wb") as op:
                pickle.dump(self.meta_data, op, 2)

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
                        "md5": md5[1],
                        "md5_path": fn,
                        "scb_comment": (
                            entry.scb_comment if hasattr(entry, "scb_comment") else ""
                        ),
                        "id": (
                            entry.scb_id if hasattr(entry, "scb_id") else None
                        ),  # this allows you to permanently tie it to an entry in SCBs databases
                        "browser_options": entry.gbrowse_options,
                        "has_splicing": entry.has_spliced_reads,
                        "genome_label": extract_genome_label(entry),
                    }
                )
                self.deps.append(md5[0])


def get_md5(input_name, input_job):
    input_name = Path(input_name)
    target_name = input_name.with_name(input_name.name + ".md5sum")

    def gen(target_filename):
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
