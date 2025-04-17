from fastapi import FastAPI, File, UploadFile, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
import os
import subprocess
import shutil
import time
import csv

app = FastAPI()
os.makedirs("uploads", exist_ok=True)
app.mount("/static", StaticFiles(directory="."), name="static")


@app.post("/download")
async def download_data(data: dict):
    accession = data.get("accession")
    if not accession:
        raise HTTPException(status_code=400, detail="未提供 NCBI 编号")
    output_zip = f"uploads/{accession}.zip"

    try:
        subprocess.run(
            ["datasets", "download", "genome", "accession", accession, "--filename", output_zip],
            check=True
        )
        subprocess.run(["unzip", "-o", output_zip, "-d", f"uploads/{accession}"], check=True)
        for root, _, files in os.walk(f"uploads/{accession}"):
            for file in files:
                if file.endswith(".fna"):
                    return {"filePath": os.path.join(root, file)}
        raise HTTPException(status_code=400, detail="未找到 FASTA 文件")
    except subprocess.CalledProcessError:
        raise HTTPException(status_code=500, detail="下载失败")


@app.post("/upload")
async def upload_file(files: list[UploadFile] = File(...)):
    file_paths = []
    for file in files:
        if not (file.filename.endswith(".fasta") or file.filename.endswith(".fna")):
            raise HTTPException(status_code=400, detail="仅支持 .fasta 或 .fna 文件")
        base, ext = os.path.splitext(file.filename)
        counter = 1
        file_path = f"uploads/{file.filename}"
        while os.path.exists(file_path):
            file_path = f"uploads/{base}_{counter}{ext}"
            counter += 1
        with open(file_path, "wb") as f:
            shutil.copyfileobj(file.file, f)
        with open(file_path, "r") as f:
            if not f.readline().strip().startswith(">"):
                os.remove(file_path)
                raise HTTPException(status_code=400, detail="无效的 FASTA 文件")
        file_paths.append(file_path)
    return {"filePaths": file_paths}


@app.get("/analyze/virulence")
async def analyze_virulence(file: str):
    if not os.path.exists(file):
        raise HTTPException(status_code=400, detail=f"文件不存在: {file}")

    try:
        result = subprocess.run(
            [
                "docker", "run", "--rm",
                "-v", f"{os.getcwd()}:/data",
                "staphb/abricate",
                "abricate", "--db", "vfdb", f"/data/{file}"
            ],
            capture_output=True,
            text=True,
            check=True
        )
        return {"virulence": result.stdout}

    except subprocess.CalledProcessError as e:
        raise HTTPException(status_code=500, detail=f"Docker 中 abricate 执行失败: {e.stderr or str(e)}")


@app.get("/analyze/resistance")
async def analyze_resistance(file: str):
    if not os.path.exists(file):
        raise HTTPException(status_code=400, detail=f"文件不存在: {file}")

    try:
        result = subprocess.run(
            [
                "docker", "run", "--rm",
                "-v", f"{os.getcwd()}:/data",
                "ncbi/amr",
                "amrfinder", "-n", f"/data/{file}"
            ],
            capture_output=True,
            text=True,
            check=True
        )
        return {"resistance": result.stdout}

    except subprocess.CalledProcessError as e:
        print("❌ AMRFinder 执行失败:", e.stderr or e.stdout)
        raise HTTPException(status_code=500, detail=f"Docker 中 AMRFinder 执行失败: {e.stderr or str(e)}")

@app.post("/analyze/nucmer")
async def analyze_nucmer(data: dict):
    files = data.get("files", [])
    if len(files) != 2:
        raise HTTPException(status_code=400, detail="nucmer 比对需正好2个文件")

    try:
        file1, file2 = files
        prefix = "uploads/nucmer_result"
        delta_file = f"{prefix}.delta"
        coords_file = f"{prefix}.coords"

        subprocess.run(["nucmer", "--maxmatch", file1, file2, "-p", prefix], check=True)
        subprocess.run(["delta-filter", "-1", delta_file], stdout=open(f"{prefix}.filtered.delta", "w"), check=True)
        subprocess.run(["show-coords", "-rcl", f"{prefix}.filtered.delta"], stdout=open(coords_file, "w"), check=True)

        with open(coords_file, "r") as f:
            coords_result = f.read()

        return {"nucmer": coords_result}

    except subprocess.CalledProcessError as e:
        raise HTTPException(status_code=500, detail=f"Nucmer 比对失败: {e.stderr or str(e)}")



@app.post("/analyze/alignment")
async def analyze_alignment(data: dict):
    files = data.get("files", [])
    if len(files) != 2:
        raise HTTPException(status_code=400, detail="全基因组比对需正好2个文件")
    try:
        for file in files:
            if not os.path.exists(file):
                raise HTTPException(status_code=400, detail=f"文件不存在: {file}")

        result = subprocess.run(
            ["minimap2", "-x", "asm5", files[0], files[1]],
            capture_output=True,
            text=True,
            check=True
        )

        # 保存 PAF 文件
        with open("uploads/alignment.paf", "w") as f:
            f.write(result.stdout)

        # 同时生成 CSV 文件（不筛选）
        csv_path = "uploads/alignment.csv"
        with open(csv_path, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([
                "query_id", "query_len", "query_start", "query_end", "strand",
                "ref_id", "ref_len", "ref_start", "ref_end",
                "match_len", "block_len", "mapq"
            ])
            for line in result.stdout.strip().splitlines():
                fields = line.strip().split()
                if len(fields) >= 13:
                    writer.writerow(fields[:13])

        return {"alignment": result.stdout}

    except subprocess.CalledProcessError as e:
        raise HTTPException(status_code=500, detail=f"Minimap2 失败: {e.stderr or str(e)}")


@app.get("/download/csv")
async def download_csv():
    csv_path = "uploads/alignment.csv"
    if not os.path.exists(csv_path):
        raise HTTPException(status_code=404, detail="CSV 文件不存在")
    return FileResponse(csv_path, media_type="text/csv", filename="alignment.csv")


@app.get("/")
async def serve_index():
    return FileResponse("index.html")
