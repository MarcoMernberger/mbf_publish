""" support code for talking to the oldschool bil2 /
getting lanes into gbrowse"""

import pickle
from pathlib import Path
import pypipegraph as ppg
import os


def activate(groups_allowed_to_access, modules):
    if isinstance(groups_allowed_to_access, str):
        groups_allowed_to_access = [groups_allowed_to_access]
    if not isinstance(modules, list):
        raise ValueError("modules must be a list")
    return BilRegistration(groups_allowed_to_access, modules).jobs


class BilRegistration:
    def __init__(self, groups_allowed_to_access, modules):
        self.modules = []
        self.add_module(FolderBrowser("results", capture_changes=False))
        self.add_module(Runscript())
        for m in modules:
            self.add_module(m)
        self.jobs = self.dump()
        self.jobs.append(self.allow_access(groups_allowed_to_access))
        self.register(self.jobs)

    def add_module(self, module):
        """Add a bil module for later dumping"""
        fn = module.get_filename()
        if not fn.endswith(".module"):
            raise ValueError("BIL module filenames need to end with .module")
        if not fn.startswith("web/modules"):
            raise ValueError("BIL module filenames need be placed in web/modules")
        self.modules.append(module)

    def dump(self):
        """Dump bil 'magic' files.

        Basically, each module can dump a file that will tell the websystem
        which modules to provide for this experiment"""
        Path("web/modules").mkdir(exist_ok=True, parents=True)
        jobs = []
        dbf = self.dump_browser_file()
        for module in self.modules:
            job = ppg.FileGeneratingJob(module.get_filename(), module.dump)
            job.depends_on(module.get_dependencies())
            if hasattr(module, "get_browser_tracks"):
                dbf.depends_on(module.get_dependencies())
            jobs.append(job)
        jobs.append(dbf)
        return jobs

    def allow_access(self, groups):
        """write the permisisons file that says which groups may access
        this project"""
        filename = Path("web/permissions.dat")

        def do_dump(groups=groups):
            if not hasattr(groups, "__iter__"):
                groups = [groups]
            filename.parent.mkdir(exist_ok=True)
            filename.write_text("\n".join(groups) + "\n")

        return ppg.FileGeneratingJob(filename, do_dump).depends_on(
            ppg.ParameterInvariant(filename, groups)
        )

    def register(self, prerequisites):
        """Call the webservice and register this project"""
        sentinel_file = Path("web/registration_sentinel")

        def do_register():
            project_name = os.environ["ANYSNAKE_PROJECT_PATH"]
            project_name = project_name[project_name.rfind("/") + 1 :]
            print("registration for", project_name)
            import requests

            data = {"project_name": project_name}
            auth = requests.auth.HTTPBasicAuth("feed", "feed")
            print(
                requests.get(
                    "http://mbf.imt.uni-marburg.de/bil2/register?",
                    params=data,
                    auth=auth,
                )
            )
            print(
                requests.get(
                    "http://mbf.imt.uni-marburg.de/bil2/gbrowse_dump", auth=auth
                )
            )

            sentinel_file.write_text("Done")

        return ppg.FileGeneratingJob(sentinel_file, do_register).depends_on(
            prerequisites
        )

    def dump_browser_file(self):
        browser_path = Path("web") / "browser.tracks"

        def do_dump_browser_file():
            tracks = []
            for module in self.modules:
                print(module)
                if hasattr(module, "get_browser_tracks"):
                    tracks.extend(module.get_browser_tracks())
            with open(browser_path, "wb") as op:
                pickle.dump(tracks, op, protocol=2)

        return ppg.FileGeneratingJob(browser_path, do_dump_browser_file)


class FolderBrowser:
    def __init__(self, folder, capture_changes=False):
        self.folder = folder
        if capture_changes:
            raise ValueError("mbf_publish.bil no longer supports capture_changes")

    def get_filename(self):
        return os.path.join(
            "web",
            "modules",
            "folder_" + self.folder.replace("/", "_").replace("\\", "_") + ".module",
        )

    def dump(self):
        data = {
            "type": "folder",
            "name": "Folder %s" % self.folder,
            "folder": self.folder,
        }
        op = open(self.get_filename(), "wb")
        pickle.dump(data, op, protocol=2)
        op.close()

    def get_dependencies(self):
        return []


class Runscript:
    """Add the run.py as a sepearate module"""

    def __init__(self, run_file="run.py"):
        self.run_file = run_file

    def get_filename(self):
        return "web/modules/runscript_%s_.module" % (self.run_file,)

    def dump(self):
        op = open(self.get_filename(), "wb")
        data = {
            "type": "sourcecode",
            "filename": self.run_file,
        }
        pickle.dump(data, op, protocol=2)
        op.close()

    def get_dependencies(self):
        return []


class AlignedLanes:
    """A listing of aligned lanes"""

    def __init__(self, lanes, browser_lanes=None):
        """Pass in either a list of lanes to show (also in the browser),
        or None ,then all lanes created will be included
        """
        if isinstance(lanes, dict):
            lanes = list(lanes.values())
        self.lanes = lanes
        self.browser_lanes = browser_lanes
        for l in self._get_lanes():
            if ":" in l.name:
                raise ValueError("Genome browser does not support : in names")

    def get_filename(self):
        return "web/modules/chipseq_aligned_lanes.module"

    def _get_lanes(self):
        return self.lanes

    def _get_browser_lanes(self):
        if self.browser_lanes is not None:
            return self.browser_lanes
        else:
            return self._get_lanes()

    def dump(self):
        data = {"name": [], "bam": [], "stats": [], "reference": []}
        for lane in self._get_lanes():
            data["name"].append(lane.name)
            data["bam"].append(lane.get_bam_names()[0])
            if hasattr(lane, "alignment_report"):
                data["stats"].append(lane.alignment_report)
            else:
                d = {}
                d["unique_perfect"], d["unique_mismatched"] = (
                    lane.mapped_reads(),
                    0,
                )
                data["stats"].append(d)
            data["reference"].append(lane.genome.name)

        import json

        print(data)
        data = {"type": "Chipseq_AlignedLanes", "data": data}
        with open(self.get_filename(), "w") as op:
            json.dump(data, op)

    def get_dependencies(self):
        deps = []
        for lane in self._get_lanes():
            deps.append(lane.load())
            # if hasattr(lane, "read_distribution_dict"):
            # deps.append(lane.read_distribution_dict())
            # else:
            # deps.append(lane.count_aligned_reads())
        # for lane in self._get_browser_lanes():
        # deps.append(lane.dump_gbrowes_adjustments())
        deps.append(
            ppg.ParameterInvariant(
                self.get_filename(),
                tuple(sorted([x.name for x in self._get_lanes()]))
                + tuple(sorted([x.name for x in self._get_browser_lanes()])),
            )
        )
        return deps

    def get_browser_tracks(self):
        res = []
        for lane in self._get_browser_lanes():
            if hasattr(lane, "gbrowse_options"):
                options = lane.gbrowse_options
            else:
                options = {}
            res.append(
                ("bam", lane.genome.name, lane.name, lane.get_bam_names()[0], options,)
            )
        return res


class GenomicRegions:
    """A listing of genomic regions, with download option"""

    def __init__(self, regions):
        if isinstance(regions, dict):
            regions = list(regions.values())
        self.regions = regions
        for r in regions:
            if ":" in r.name:
                raise ValueError("Genome browser does not support : in names")

    def get_filename(self):
        return "web/modules/chipseq_peaksets.module"

    def _get_regions(self):
        return self.regions

    def to_outside_path(f, str_path):
        return str_path.replace("/project/", os.environ["ANYSNAKE_PROJECT_PATH"] + "/")

    def dump(self):
        data = {
            "name": [],
            "filename": [],
            "peak_count": [],
        }
        for peakset in self._get_regions():
            name = peakset.name
            filepath = str(peakset.write()[1])
            data["name"].append(name)
            if os.path.exists(filepath):
                data["filename"].append(self.to_outside_path(filepath))
            else:
                data["filename"].append("")
            data["peak_count"].append(len(peakset.df))
        if data:
            import json

            data = {"type": "Chipseq_Peaksets", "data": data}
            with open(self.get_filename(), "w") as op:
                json.dump(data, op)

    def get_dependencies(self):
        res = []
        for peakset in self._get_regions():
            res.append(peakset.load())
            res.append(peakset.write_bigbed()[0])
            res.append(peakset.write()[0])
        res.append(
            ppg.ParameterInvariant(
                self.get_filename(),
                tuple(sorted([x.name for x in self._get_regions()])),
            )
        )
        return res

    def get_browser_tracks(self):
        res = []
        for region in self._get_regions():
            res.append(
                (
                    "bigbed",
                    region.genome.name,
                    region.name,
                    self.to_outside_path(str(region.write_bigbed()[1])),
                    {},
                )
            )
        return res
