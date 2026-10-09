"""Parse VCV XML and index ClinVar's explicit submitted-trait mappings.

The temporary SQLite index keeps full-release memory use bounded. It is an
internal parser working file, not part of the phenotype database or JSONL.
"""

import gzip
import sqlite3
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from time import monotonic


# Shared default for VCV indexing and RCV parsing; counts input records.
PROGRESS_INTERVAL = 100_000


def format_elapsed(started):
    hours, remainder = divmod(int(monotonic() - started), 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"


def normalise_name(value):
    return " ".join(value.split()).casefold()


def mapping_key(kind, reference, value):
    if kind == "Name":
        value = normalise_name(value)
    else:
        value = value.removeprefix(reference + ":")
    return kind, reference.casefold(), value


def trait_keys(trait):
    """Keep raw trait context: exported identifiers may omit some useful XRefs."""
    result = set()
    for name in trait.findall("Name/ElementValue"):
        if name.text and name.text.strip():
            result.add(mapping_key("Name", name.get("Type", ""), name.text.strip()))
    for path in ("XRef", "Name/XRef", "Symbol/XRef", "AttributeSet/XRef"):
        for xref in trait.findall(path):
            database, identifier = xref.get("DB"), xref.get("ID")
            if database and identifier:
                result.add(mapping_key("XRef", database, identifier))
    return result


class VcvMappings:
    def __init__(self, path, *, directory):
        self.path = Path(path)
        self.directory = directory
        self.release_date = None
        self.n_records = 0
        self._temporary = None
        self._db = None

    def __enter__(self):
        self._temporary = tempfile.TemporaryDirectory(prefix="clinvar-vcv-", dir=self.directory)
        try:
            self._db = sqlite3.connect(Path(self._temporary.name) / "mappings.sqlite")
            self._db.execute("CREATE TABLE rcv_versions (vcv TEXT, rcv TEXT)")
            self._db.execute("""CREATE TABLE mappings (
                vcv TEXT, scv TEXT, trait_type TEXT, kind TEXT, reference TEXT,
                value TEXT, cui TEXT, name TEXT
            )""")
            self._build()
            return self
        except BaseException:
            self.__exit__(None, None, None)
            raise

    def __exit__(self, *args):
        if self._db is not None:
            self._db.close()
        if self._temporary is not None:
            self._temporary.cleanup()

    def _progress(self, stage, started):
        print(f"[Parsing ClinVar VCV] {stage}: {self.n_records:,} records processed; "
              f"elapsed {format_elapsed(started)}",
              file=sys.stderr, flush=True)

    def _build(self):
        started = monotonic()
        print(f"[Parsing ClinVar VCV] Starting: {self.path}", file=sys.stderr, flush=True)
        opener = gzip.open if self.path.suffix.lower() == ".gz" else open
        root = None
        with opener(self.path, "rb") as stream, self._db:
            for event, element in ET.iterparse(stream, events=("start", "end")):
                if root is None and event == "start":
                    if element.tag not in {"ClinVarVariationRelease", "ClinVarResult-Set"}:
                        raise ValueError(f"Expected VCV XML root, found {element.tag}")
                    root = element
                    self.release_date = root.get("ReleaseDate")
                    continue
                if event != "end" or element.tag != "VariationArchive":
                    continue
                self.n_records += 1
                if element.findtext("RecordStatus") in (None, "current"):
                    self._add_record(element)
                element.clear()
                root.clear()
                if self.n_records % PROGRESS_INTERVAL == 0:
                    self._progress("Progress", started)
            if root is None:
                raise ValueError("No VCV XML root found")
            self._progress("Finalising index", started)
            self._db.execute("CREATE INDEX rcv_version_key ON rcv_versions (vcv, rcv)")
            self._db.execute("CREATE INDEX mapping_key ON mappings (vcv, scv, kind, reference, value)")
        self._progress("Complete", started)

    def _add_record(self, archive):
        vcv = archive.get("Accession")
        classified = archive.find("ClassifiedRecord")
        if not vcv or classified is None:
            return
        rcv_rows = []
        for accession in classified.findall("RCVList/RCVAccession"):
            identifier, version = accession.get("Accession"), accession.get("Version")
            if identifier and version:
                rcv_rows.append((vcv, f"{identifier}.{version}"))
        self._db.executemany("INSERT INTO rcv_versions VALUES (?, ?)", rcv_rows)

        # Join using ClinicalAssertion/@ID, not the digits in an SCV accession.
        # The internal assertion ID can differ from the accession's number.
        scvs = {}
        for assertion in classified.findall("ClinicalAssertionList/ClinicalAssertion"):
            if assertion.findtext("RecordStatus") not in (None, "current"):
                continue
            accession = assertion.find("ClinVarAccession")
            if (accession is not None and accession.get("Accession")
                    and accession.get("Version") and assertion.get("ID")):
                identifier, version = accession.get("Accession"), accession.get("Version")
                scvs[assertion.get("ID")] = f"{identifier}.{version}"
        rows = []
        for mapping in classified.findall("TraitMappingList/TraitMapping"):
            assertion_id = mapping.get("ClinicalAssertionID")
            scv = scvs.get(assertion_id)
            medgen = mapping.find("MedGen")
            kind, reference, value = (mapping.get(attr, "") for attr in
                                      ("MappingType", "MappingRef", "MappingValue"))
            if not scv or medgen is None or kind not in {"Name", "XRef"} or not value:
                continue
            kind, reference, value = mapping_key(kind, reference, value)
            rows.append((vcv, scv, mapping.get("TraitType", ""), kind, reference, value,
                         medgen.get("CUI"), medgen.get("Name")))
        self._db.executemany("INSERT INTO mappings VALUES (?, ?, ?, ?, ?, ?, ?, ?)", rows)

    def has_rcv_version(self, vcv, rcv):
        if not vcv:
            return False
        return self._db.execute("SELECT 1 FROM rcv_versions WHERE vcv = ? AND rcv = ? LIMIT 1",
                                (vcv.split(".")[0], rcv)).fetchone() is not None

    def find(self, vcv, scv, trait):
        matches = set()
        for kind, reference, value in trait_keys(trait):
            rows = self._db.execute("""SELECT trait_type, cui, name
                FROM mappings WHERE vcv = ? AND scv = ? AND kind = ? AND reference = ? AND value = ?""",
                (vcv.split(".")[0], scv, kind, reference, value))
            for trait_type, cui, name in rows:
                if trait_type and trait.get("Type") and trait_type != trait.get("Type"):
                    continue
                matches.add((cui, name))
        return sorted(matches, key=lambda row: tuple(value or "" for value in row))
