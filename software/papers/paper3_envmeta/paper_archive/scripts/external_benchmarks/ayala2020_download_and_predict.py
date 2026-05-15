"""
Ayala-Muñoz 2020 Microorganisms 13 MAGs 下载 + ORF 预测 一键脚本。

用途: 对照实验 Arm C2 — 准备 GhostKOALA 上传文件。

依赖（用户首次跑前安装一次）:
    pip install pyrodigal requests

用法:
    python tools/external_benchmarks/ayala2020_download_and_predict.py \\
        --out paper/benchmarks/external/ayala_2020_pitlake/input_data_local

输出:
    input_data_local/
    ├── mag_fasta/                # 13 个解压后的 MAG 核酸 fasta
    │   ├── A_NAN_12.fna
    │   └── ...
    ├── mag_proteins/             # 13 个 prodigal 预测 蛋白 fasta
    │   ├── A_NAN_12.faa
    │   └── ...
    └── all_mags_proteins.faa     # 合并 → 上传 GhostKOALA（< 200 MB）

之后用户操作:
    1. 打开 https://www.kegg.jp/ghostkoala/
    2. 上传 all_mags_proteins.faa（输入邮箱）
    3. 选 "Genus_prokaryotes + family_eukaryotes" KEGG GENES set
    4. 等 4-24 h email 通知（异步）
    5. 下载结果 user.out → 放到 input_data_local/ghostkoala_user.out
    6. 跑 ayala2020_reshape.py（待 KO 拿到后写）
"""

from __future__ import annotations

import argparse
import gzip
import io
import sys
import time
from pathlib import Path
from urllib.request import urlopen, Request

# 13 MAGs 来自 PRJNA646106（NCBI BioProject 截图确认）
# Format: (GCA_accession, isolate_name, taxonomy)
MAGS = [
    ("GCA_014874275.1", "B_ACI_09", "Acidobacteriota_bacterium"),
    ("GCA_014874305.1", "B_PAT_13", "Colwellibacteriota_bacterium"),
    ("GCA_014874255.1", "A_NAN_12", "Woesearchaeota_archaeon"),
    ("GCA_014874385.1", "B_CHL_03", "Dehalococcoidia_bacterium"),
    ("GCA_014874295.1", "A_MIC_10", "Microcaldota_archaeon"),
    ("GCA_014874415.1", "A_CRE_07", "Nitrososphaerales_archaeon"),
    ("GCA_014874185.1", "B_NIT_04", "Nitrospirota_bacterium"),
    ("GCA_014874425.1", "B_PRO_05", "Syntrophaceae_bacterium"),
    ("GCA_014874265.1", "B_ACT_11", "Thermoleophilia_bacterium"),
    ("GCA_014874175.1", "B_ACT_02", "Thermoleophilia_bacterium"),
    ("GCA_014874235.1", "A_EUR_01", "Thermoplasmatales_archaeon"),
    ("GCA_014874365.1", "A_EUR_06", "Thermoplasmatales_archaeon"),
    ("GCA_014874355.1", "B_DOR_08", "bacterium"),
]


def gca_to_ftp_url(gca: str) -> str:
    """Build NCBI FTP URL for a GCA accession.

    GCA_014874275.1 → https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/874/275/
    然后需要列目录找到具体 _genomic.fna.gz 文件名。
    """
    parts = gca.replace("GCA_", "").split(".")[0]
    a, b, c = parts[0:3], parts[3:6], parts[6:9]
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{a}/{b}/{c}/"


def find_genome_file(base_url: str, gca: str) -> str | None:
    """List FTP dir, find {accession}_*/...genomic.fna.gz file URL."""
    try:
        req = Request(base_url, headers={"User-Agent": "Mozilla/5.0"})
        with urlopen(req, timeout=30) as r:
            content = r.read().decode("utf-8", errors="replace")
        # Find subdirectory link starting with GCA accession
        import re
        m = re.search(rf'href="({re.escape(gca)}[^"]*)/"', content)
        if not m:
            return None
        sub = m.group(1)
        # Build URL to specific file
        return f"{base_url}{sub}/{sub}_genomic.fna.gz"
    except Exception as e:
        print(f"    [ERROR] list ftp dir failed: {e}")
        return None


def download_mag(gca: str, name: str, out_dir: Path) -> Path | None:
    """Download genomic.fna.gz for one MAG, decompress to out_dir/{name}.fna"""
    out_path = out_dir / f"{name}.fna"
    if out_path.exists():
        print(f"  [SKIP] {name}.fna already exists ({out_path.stat().st_size//1024} KB)")
        return out_path

    base_url = gca_to_ftp_url(gca)
    file_url = find_genome_file(base_url, gca)
    if not file_url:
        print(f"  [ERROR] Could not find genomic.fna.gz for {gca} at {base_url}")
        return None

    try:
        print(f"  Downloading {name} <- {file_url.split('/')[-1]}")
        req = Request(file_url, headers={"User-Agent": "Mozilla/5.0"})
        with urlopen(req, timeout=120) as r:
            gz_data = r.read()
        # Decompress
        with gzip.GzipFile(fileobj=io.BytesIO(gz_data)) as gz:
            content = gz.read()
        out_path.write_bytes(content)
        print(f"    → {out_path.name} ({len(content)//1024} KB)")
        return out_path
    except Exception as e:
        print(f"  [ERROR] download {gca} failed: {e}")
        return None


def predict_orfs(fna_path: Path, faa_path: Path, name: str) -> int:
    """Predict ORFs using pyrodigal, write protein fasta."""
    try:
        import pyrodigal
    except ImportError:
        print("ERROR: pyrodigal not installed. Run: pip install pyrodigal")
        sys.exit(1)

    # MAG 通常 single-genome mode, 但也有 metagenome 模式 — 我们用 meta 更鲁棒
    orf_finder = pyrodigal.GeneFinder(meta=True)

    n_proteins = 0
    sequences = []
    # Parse fna
    with open(fna_path) as f:
        contig_id = None
        seq_buf: list[str] = []
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if contig_id and seq_buf:
                    sequences.append((contig_id, "".join(seq_buf)))
                contig_id = line[1:].split()[0]
                seq_buf = []
            else:
                seq_buf.append(line)
        if contig_id and seq_buf:
            sequences.append((contig_id, "".join(seq_buf)))

    with open(faa_path, "w") as out:
        for contig_id, seq in sequences:
            try:
                genes = orf_finder.find_genes(seq.encode("ascii", errors="ignore"))
            except Exception as e:
                print(f"    [WARN] {name}/{contig_id} pyrodigal failed: {e}")
                continue
            for i, gene in enumerate(genes, 1):
                prot = gene.translate()
                # Header format: {MAG_name}|{contig_id}|gene_{i}
                out.write(f">{name}|{contig_id}|gene_{i}\n{prot}\n")
                n_proteins += 1
    return n_proteins


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True, help="output directory for input_data_local")
    args = p.parse_args()

    out_root = Path(args.out)
    fna_dir = out_root / "mag_fasta"
    faa_dir = out_root / "mag_proteins"
    fna_dir.mkdir(parents=True, exist_ok=True)
    faa_dir.mkdir(parents=True, exist_ok=True)

    print(f"Step 1/3: Download 13 MAG fasta from NCBI FTP")
    print(f"  → {fna_dir}\n")

    downloaded: list[tuple[str, Path]] = []
    for gca, name, tax in MAGS:
        path = download_mag(gca, name, fna_dir)
        if path:
            downloaded.append((name, path))
        time.sleep(0.5)  # be polite to NCBI FTP
    print(f"\n  {len(downloaded)} / {len(MAGS)} MAGs downloaded\n")

    if not downloaded:
        print("FATAL: no MAGs downloaded; aborting")
        sys.exit(1)

    print(f"Step 2/3: Predict ORFs with pyrodigal")
    print(f"  → {faa_dir}\n")
    total_proteins = 0
    for name, fna_path in downloaded:
        faa_path = faa_dir / f"{name}.faa"
        if faa_path.exists():
            print(f"  [SKIP] {name}.faa exists")
            with open(faa_path) as f:
                n = sum(1 for line in f if line.startswith(">"))
            total_proteins += n
            continue
        print(f"  Predicting ORFs for {name}...")
        n = predict_orfs(fna_path, faa_path, name)
        print(f"    → {n} proteins")
        total_proteins += n
    print(f"\n  Total proteins: {total_proteins}\n")

    print(f"Step 3/3: Concatenate all proteins for GhostKOALA upload")
    combined = out_root / "all_mags_proteins.faa"
    with open(combined, "w") as out:
        for name, _ in downloaded:
            faa_path = faa_dir / f"{name}.faa"
            if faa_path.exists():
                out.write(faa_path.read_text())
    size_mb = combined.stat().st_size / (1024 * 1024)
    print(f"  → {combined.name} ({size_mb:.1f} MB)")
    if size_mb > 300:
        print(f"  ⚠️ WARNING: file > 300 MB GhostKOALA limit. Split into batches.")
    else:
        print(f"  ✅ < 300 MB (GhostKOALA accepts up to 300 MB single submission)")

    print(f"\nDone. Next steps for user:")
    print(f"  1. Open https://www.kegg.jp/ghostkoala/")
    print(f"  2. Upload {combined}")
    print(f"  3. Select 'Genus_prokaryotes + family_eukaryotes' (KEGG GENES set)")
    print(f"  4. Wait 4-24 h for email with results download link")
    print(f"  5. Download user.out → save to:")
    print(f"     {out_root / 'ghostkoala_user.out'}")
    print(f"  6. Then run: python tools/external_benchmarks/ayala2020_reshape.py")


if __name__ == "__main__":
    main()
