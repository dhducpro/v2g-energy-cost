"""Recreate the paper's study and July instances from the pinned processed Caltech export."""

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timedelta, timezone
import hashlib
import io
import json
import math
from pathlib import Path
from zipfile import ZipFile
from zoneinfo import ZoneInfo


LOCAL = ZoneInfo("America/Los_Angeles")
UTC = timezone.utc
ARCHIVE_SHA256 = "484b607820ca818757eef6b9ecc7c09b8d68149a0cb4123e66100ae23f7c0a3c"
INPUT_PREFIX = "replication/cohort_reconstruction/inputs/"
INPUT_SHA256 = {
    "caltech.csv": "c7ba7126156ffb99a926258f00abe0c15d5a33341c4580a042189d04b600b861",
    "frozen_selected_sessions.csv": "c769b8298a7ed11e5a25168cde93f69830ed67b03d2b0e2866a22ac4c17f2504",
    "original_excluded_sessions.csv": "7a38d1bd0f5b8c2105279c450282237bccbf202c590de3219486177a98ce90c5",
}


def read_tables(archive):
    if hashlib.sha256(archive.read_bytes()).hexdigest() != ARCHIVE_SHA256:
        raise ValueError("Archive does not match the published input commit.")
    with ZipFile(archive) as source:
        tables = {}
        for filename, checksum in INPUT_SHA256.items():
            raw = source.read(INPUT_PREFIX + filename)
            if hashlib.sha256(raw).hexdigest() != checksum:
                raise ValueError(f"Unexpected source contents in {filename}")
            tables[filename] = list(csv.DictReader(io.StringIO(raw.decode("utf-8"))))
    return tables


def local_time(value):
    moment = datetime.fromisoformat(value.replace("Z", "+00:00"))
    return moment.replace(tzinfo=LOCAL) if moment.tzinfo is None else moment.astimezone(LOCAL)


def hours(start, end):
    return (end.astimezone(UTC) - start.astimezone(UTC)).total_seconds() / 3600


def key_for(start, end, energy):
    return "|".join((start.astimezone(UTC).isoformat(), end.astimezone(UTC).isoformat(),
                     format(energy, ".12g")))


def digest(token):
    return hashlib.sha256(("20261114|" + token).encode("utf-8")).hexdigest()


def pool_from_tables(tables):
    excluded = {r["physical_key"] for r in tables["frozen_selected_sessions.csv"]}
    excluded.update(key_for(local_time(r["connection_local"]), local_time(r["disconnect_local"]),
                            float(r["reported_delivered_kwh_as_required_battery_increment"]))
                    for r in tables["original_excluded_sessions.csv"])
    if len(excluded) != 162 or len(tables["caltech.csv"]) != 27558:
        raise ValueError("Historical source record counts changed")
    pool, july_screen, seen = defaultdict(list), Counter(), set()
    for row in tables["caltech.csv"]:
        try:
            start, end = local_time(row["Start"]), local_time(row["End"])
            energy = float(row["Energy"])
        except (ValueError, TypeError):
            continue
        day = start.date().isoformat()
        if not "2019-04-01" <= day < "2019-10-01":
            continue
        key = key_for(start, end, energy)
        def assign(field, options):
            token = f"observed|{field}|{key}"
            return options[int(digest(token), 16) % len(options)]
        origin = start.replace(hour=0, minute=0, second=0, microsecond=0)
        last = origin + timedelta(days=1)
        pc = assign("Pc", (3.6, 7.2, 11.0))
        if not math.isfinite(energy) or energy < 2:
            status = "energy_below_screen"
        elif end <= start:
            status = "nonpositive_duration"
        elif end > last:
            status = "cross_day"
        elif hours(origin, last) != 24:
            status = "non_24hour_day"
        elif key in seen:
            status = "duplicate"
        elif key in excluded:
            status = "prior_fleet_overlap"
        elif energy > pc * hours(start, end) + 1e-9:
            status = "charger_deliverability"
        else:
            status = "eligible"
        seen.add(key)
        if day.startswith("2019-07-"):
            july_screen[status] += 1
        if status != "eligible":
            continue
        initial = assign("E0", (15., 25., 35.))
        eta_c = assign("eta_c", (.90, .93, .96))
        target = initial + eta_c * energy
        pool[day].append(dict(
            key=key, arrival=hours(origin, start) * 3600,
            departure=hours(origin, end) * 3600, E0=initial, Emin=0.,
            Etarget=target, Emax=target + assign("headroom", (0., .5, 5., 15.)),
            Pc=pc, Pd=assign("Pd", (2., 4., 7.2)), eta_c=eta_c,
            eta_d=assign("eta_d", (.88, .92, .96))))
    return pool, july_screen


def fleet(name, day, sessions):
    vehicles = sorted(sessions, key=lambda v: (v["arrival"], v["departure"], v["key"]))
    edges = sorted({float(t * 3600) for t in range(25)}
                   | {v[k] for v in vehicles for k in ("arrival", "departure")})
    intervals = list(zip(edges, edges[1:]))
    nt = len(intervals)
    model = {
        "name": name, "N": len(vehicles), "T": nt,
        "availability": [[int(v["arrival"] <= (lo+hi)/2 < v["departure"])
                          for lo, hi in intervals] for v in vehicles],
        "duration": [(hi-lo)/3600 for lo, hi in intervals],
        "boundaries_seconds": edges, "L": [10.]*nt, "Ulo": [1.]*nt, "Uhi": [1.]*nt,
    }
    for field in ("E0", "Emin", "Emax", "Etarget", "Pc", "Pd", "eta_c", "eta_d"):
        model[field] = [v[field] for v in vehicles]
    model["metadata"] = {"horizon_date": day}
    return model


def build_cases(pool):
    cases = {}
    sizes = (6, 12, 24, 6, 12, 24)
    for month in range(4, 10):
        idx = 0
        days = sorted(d for d in pool if d.startswith(f"2019-{month:02d}-"))
        for day in days:
            n = sizes[idx]
            if len(pool[day]) < n:
                continue
            ranked = sorted(pool[day], key=lambda v: (digest(f"observed|selection|{day}|{v['key']}"),
                                                      v["key"]))
            name = f"A2019{month:02d}_{idx+1:02d}_N{n:02d}"
            cases[f"inputs/study/{name}.json"] = fleet(name, day, ranked[:n])
            idx += 1
            if idx == len(sizes):
                break
        if idx != len(sizes):
            raise ValueError(f"Insufficient eligible dates in month {month}")
    for i in range(1, 32):
        day = f"2019-07-{i:02d}"
        name = "J" + day.replace("-", "")
        cases[f"inputs/july/{name}.json"] = fleet(name, day, pool[day])
    return cases


def main():
    root = Path(__file__).resolve().parent.parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--archive", type=Path, required=True,
                    help="Pinned study-data.zip from the documented repository commit")
    ap.add_argument("--manifest", type=Path, default=root / "inputs_sha256.csv")
    ap.add_argument("--out", type=Path,
                    help="Optional new directory for verified files; existing files are preserved")
    args = ap.parse_args()
    with args.manifest.open(newline="", encoding="utf-8") as f:
        pinned = {r["relative_path"]: r["sha256"] for r in csv.DictReader(f)}
    pool, july_screen = pool_from_tables(read_tables(args.archive))
    if july_screen != {"eligible": 578, "energy_below_screen": 158,
                       "cross_day": 39, "charger_deliverability": 36}:
        raise ValueError("July screening counts differ from the study")
    cases = build_cases(pool)
    if len(cases) != 67 or set(cases) != set(pinned):
        raise ValueError("Generated file names do not match the 67 input files in the manifest")
    checked = {}
    for name, model in cases.items():
        # Store exactly the canonical byte representation of the numerical archive.
        serialized = (json.dumps(model, separators=(",", ":"), allow_nan=False) + "\r\n").encode("utf-8")
        if hashlib.sha256(serialized).hexdigest() != pinned[name]:
            raise ValueError(f"Generated instance differs from frozen input: {name}")
        checked[name] = serialized
    if args.out is not None:
        if any((args.out / name).exists() for name in checked):
            raise FileExistsError("Choose a new output folder; an output file already exists.")
        for name, data in checked.items():
            target = args.out / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
    for kind in ("study", "july"):
        members = [m for name,m in cases.items() if name.startswith(f"inputs/{kind}/")]
        print(f"{kind}: {len(members)} fleet-days, {sum(m['N'] for m in members)} selected sessions")
    print(f"July eligibility: {sum(july_screen.values())} arrivals; {dict(sorted(july_screen.items()))}")
    print(f"Exact verification: {len(checked)}/{len(pinned)} archived SHA-256 digests matched.")


if __name__ == "__main__":
    main()
