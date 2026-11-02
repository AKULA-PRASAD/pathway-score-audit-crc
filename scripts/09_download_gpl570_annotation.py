import urllib.request
from pathlib import Path

url = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    "?targ=self&acc=GPL570&form=text&view=full"
)

output_dir = Path("data/raw/gpl570")
output_dir.mkdir(parents=True, exist_ok=True)

output_file = output_dir / "GPL570_platform.soft"

request = urllib.request.Request(
    url,
    headers={
        "User-Agent": "Mozilla/5.0",
        "Accept": "text/plain"
    }
)

print("Retrieving GPL570 platform record from NCBI GEO...")

with urllib.request.urlopen(request, timeout=120) as response:
    content = response.read()

if not content:
    raise ValueError("NCBI returned an empty response.")

text_content = content.decode("utf-8", errors="replace")

if "^PLATFORM = GPL570" not in text_content:
    raise ValueError(
        "The response does not contain the expected GPL570 platform record. "
        "The downloaded response will not be saved as valid annotation."
    )

output_file.write_text(text_content, encoding="utf-8")

print("Download successful.")
print("Saved to:", output_file)
print("File size (MB):", round(len(content) / (1024 * 1024), 2))
print("GPL570 platform record confirmed.")
