#!/usr/bin/env python3
import re
from pysradb.sraweb import SRAweb

GSM = [
"GSM3211196","GSM3211194","GSM3211169","GSM3211171","GSM3211174","GSM3211175",
"GSM3211177","GSM3211178","GSM3211186","GSM3211185","GSM3211180","GSM3211179",
"GSM3211183","GSM3211182","GSM3211189","GSM3211190","GSM3211193","GSM3211191"
]

db = SRAweb()
df = db.metadata("SRP151065", detailed=True)

if "library_name" not in df.columns:
    raise SystemExit("pysradb metadata does not contain library_name")

text_cols = [c for c in ["experiment_title", "sample_title", "experiment_desc"] if c in df.columns]
if not text_cols:
    raise SystemExit("pysradb metadata has no searchable sample/experiment title columns")

mapping = {}
for _, row in df.iterrows():
    text = " ".join(str(row.get(c, "")) for c in text_cols)
    hits = re.findall(r"GSM\d+", text)
    for gsm in hits:
        if gsm in GSM and gsm not in mapping:
            mapping[gsm] = str(row["library_name"])

missing = [g for g in GSM if g not in mapping]
if missing:
    raise SystemExit("Missing GSM mappings: " + ",".join(missing))

import pandas as pd
out = pd.DataFrame({"experiment_geo_accession": GSM, "library_name": [mapping[g] for g in GSM]})
out.to_csv("results/library_mapping.tsv", sep="\t", index=False)
