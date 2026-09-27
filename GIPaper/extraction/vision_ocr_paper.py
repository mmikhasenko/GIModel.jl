#!/usr/bin/env python3
"""Vision-OCR a two-column PDF into equation-preserving Markdown.

This script renders each PDF page, creates overlapping left/right column crops,
then sends a low-detail full page plus high-detail column crops to a vision
model. Outputs are cached per page so interrupted runs can resume.
"""

from __future__ import annotations

import argparse
import base64
import json
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


# The scanned paper is not distributed with the packages. Point GI_PAPER_DIR at
# a local folder holding Godfrey-Isgur-1985.pdf; OCR output is written beside it.
PAPER_DIR = Path(os.environ.get("GI_PAPER_DIR", "paper")).resolve()
PDF = PAPER_DIR / "Godfrey-Isgur-1985.pdf"
OUT_DIR = PAPER_DIR / "vision_ocr"
MODEL = "gpt-4.1"


PROMPT = """You are transcribing a scanned two-column physics paper page.

Inputs:
- Image 1 is the full page at low detail, used only for page layout.
- Image 2 is the left column crop at high detail.
- Image 3 is the right column crop at high detail.

Return clean Markdown for this page only. Do not wrap the response in a
Markdown code fence.

Rules:
- Preserve reading order: page header if meaningful, then left column, then right column.
- Preserve equations as LaTeX display math. Keep equation numbers like \\tag{A10}.
- For inline math, use LaTeX between `$...$`.
- Do not invent missing symbols. If a symbol is genuinely unclear, write `[unclear: ...]`.
- Transcribe symbol-by-symbol from the images; do not correct formulas from
  memory or from later context.
- Preserve subscripts, superscripts, bars, tildes, hats, Greek letters, primes, dots, vectors, fractions, and signs carefully.
- Preserve figure/table captions and table text if visible, but do not try to make perfect CSV tables.
- Remove obvious OCR artifacts and line-break hyphenation in prose.
- If text is cut by a crop boundary, use the full-page image to resolve it when possible.
- Do not summarize or explain the page.
"""


def run(cmd: list[str]) -> str:
    proc = subprocess.run(cmd, check=True, text=True, capture_output=True)
    return proc.stdout


def page_count(pdf: Path) -> int:
    info = run(["pdfinfo", str(pdf)])
    match = re.search(r"^Pages:\s+(\d+)", info, re.MULTILINE)
    if not match:
        raise RuntimeError("could not determine PDF page count from pdfinfo")
    return int(match.group(1))


def image_size(path: Path) -> tuple[int, int]:
    out = run(["magick", "identify", "-format", "%w %h", str(path)])
    w, h = out.split()
    return int(w), int(h)


def ensure_rendered(pdf: Path, page: int, pages_dir: Path, dpi: int) -> Path:
    pages_dir.mkdir(parents=True, exist_ok=True)
    page_path = pages_dir / f"page-{page:03d}.png"
    if page_path.exists():
        return page_path

    prefix = pages_dir / f"render-{page:03d}"
    run([
        "pdftoppm",
        "-f",
        str(page),
        "-l",
        str(page),
        "-r",
        str(dpi),
        "-png",
        str(pdf),
        str(prefix),
    ])
    rendered = pages_dir / f"render-{page:03d}-{page}.png"
    if not rendered.exists():
        matches = sorted(pages_dir.glob(f"render-{page:03d}-*.png"))
        if not matches:
            raise RuntimeError(f"pdftoppm did not create a page image for page {page}")
        rendered = matches[0]
    rendered.rename(page_path)
    return page_path


def ensure_crops(page_image: Path, crops_dir: Path, page: int) -> tuple[Path, Path]:
    crops_dir.mkdir(parents=True, exist_ok=True)
    left = crops_dir / f"page-{page:03d}-left.png"
    right = crops_dir / f"page-{page:03d}-right.png"
    if left.exists() and right.exists():
        return left, right

    w, h = image_size(page_image)
    overlap = max(80, w // 30)
    half = w // 2
    left_w = min(w, half + overlap)
    right_x = max(0, half - overlap)
    right_w = w - right_x

    run(["magick", str(page_image), "-crop", f"{left_w}x{h}+0+0", "+repage", str(left)])
    run(["magick", str(page_image), "-crop", f"{right_w}x{h}+{right_x}+0", "+repage", str(right)])
    return left, right


def data_url(path: Path) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def extract_output_text(payload: dict) -> str:
    parts: list[str] = []
    for item in payload.get("output", []):
        if item.get("type") == "message":
            for content in item.get("content", []):
                if content.get("type") == "output_text":
                    parts.append(content.get("text", ""))
    if parts:
        return "\n".join(parts).strip()
    text = payload.get("output_text")
    if isinstance(text, str):
        return text.strip()
    return ""


def clean_model_markdown(text: str) -> str:
    text = text.strip()
    fence = re.fullmatch(r"```(?:markdown|md)?\s*(.*?)\s*```", text, flags=re.DOTALL)
    if fence:
        return fence.group(1).strip()
    return text


def call_openai(
    api_key: str,
    model: str,
    page: int,
    total_pages: int,
    full_page: Path,
    left: Path,
    right: Path,
    max_output_tokens: int,
) -> tuple[str, dict]:
    body = {
        "model": model,
        "input": [
            {
                "role": "user",
                "content": [
                    {
                        "type": "input_text",
                        "text": f"{PROMPT}\nPDF page {page} of {total_pages}.",
                    },
                    {
                        "type": "input_image",
                        "image_url": data_url(full_page),
                        "detail": "low",
                    },
                    {
                        "type": "input_image",
                        "image_url": data_url(left),
                        "detail": "high",
                    },
                    {
                        "type": "input_image",
                        "image_url": data_url(right),
                        "detail": "high",
                    },
                ],
            }
        ],
        "max_output_tokens": max_output_tokens,
    }
    req = urllib.request.Request(
        "https://api.openai.com/v1/responses",
        data=json.dumps(body).encode("utf-8"),
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
        },
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=180) as response:
        payload = json.loads(response.read().decode("utf-8"))
    text = extract_output_text(payload)
    if not text:
        raise RuntimeError(f"no output text returned for page {page}")
    return clean_model_markdown(text), payload.get("usage", {})


def append_jsonl(path: Path, row: dict) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdf", type=Path, default=PDF)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--model", default=MODEL)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--start-page", type=int, default=1)
    parser.add_argument("--end-page", type=int)
    parser.add_argument("--max-output-tokens", type=int, default=4096)
    parser.add_argument("--sleep", type=float, default=0.2)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    args.pdf = args.pdf.resolve()
    args.out_dir = args.out_dir.resolve()

    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise SystemExit("OPENAI_API_KEY is not set")

    total_pages = page_count(args.pdf)
    end_page = args.end_page or total_pages
    if args.start_page < 1 or end_page > total_pages or args.start_page > end_page:
        raise SystemExit(f"invalid page range {args.start_page}-{end_page}; PDF has {total_pages} pages")

    pages_dir = args.out_dir / "page_images"
    crops_dir = args.out_dir / "column_crops"
    page_md_dir = args.out_dir / "pages"
    page_md_dir.mkdir(parents=True, exist_ok=True)
    usage_path = args.out_dir / "usage.jsonl"

    completed: list[Path] = []
    for page in range(args.start_page, end_page + 1):
        md_path = page_md_dir / f"page-{page:03d}.md"
        if md_path.exists() and not args.force:
            print(f"page {page}: cached {md_path}")
            completed.append(md_path)
            continue

        print(f"page {page}: render/crop", flush=True)
        full_page = ensure_rendered(args.pdf, page, pages_dir, args.dpi)
        left, right = ensure_crops(full_page, crops_dir, page)

        print(f"page {page}: OpenAI {args.model}", flush=True)
        try:
            text, usage = call_openai(
                api_key,
                args.model,
                page,
                total_pages,
                full_page,
                left,
                right,
                args.max_output_tokens,
            )
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"OpenAI API error on page {page}: {exc.code} {detail}") from exc

        md_path.write_text(f"<!-- PDF page {page}; model {args.model} -->\n\n{text}\n", encoding="utf-8")
        append_jsonl(
            usage_path,
            {
                "page": page,
                "model": args.model,
                "usage": usage,
                "markdown": str(md_path.relative_to(PAPER_DIR)),
                "full_page": str(full_page.relative_to(PAPER_DIR)),
                "left_crop": str(left.relative_to(PAPER_DIR)),
                "right_crop": str(right.relative_to(PAPER_DIR)),
            },
        )
        completed.append(md_path)
        if args.sleep:
            time.sleep(args.sleep)

    aggregate = args.out_dir / "godfrey_isgur_1985_vision_ocr.md"
    chunks = [
        "# Godfrey-Isgur 1985 Vision OCR\n\n",
        "Generated from rendered page/column images with a vision model. ",
        "Use the PDF and saved crops as the authority for final equation audits.\n\n",
    ]
    for page in range(1, total_pages + 1):
        md_path = page_md_dir / f"page-{page:03d}.md"
        if md_path.exists():
            chunks.append(f"## PDF Page {page}\n\n")
            body = md_path.read_text(encoding="utf-8").strip()
            body = re.sub(r"^<!--.*?-->\s*", "", body, flags=re.DOTALL)
            chunks.append(body)
            chunks.append("\n\n")
    aggregate.write_text("".join(chunks), encoding="utf-8")
    print(f"wrote {aggregate}")
    print(f"usage log: {usage_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
