"""Ensembl REST API client with retry logic and batch sequence fetching."""

import sys
import time
from typing import Any, Dict, List, Optional

import requests


BASE_URL = "https://rest.ensembl.org"
BATCH_SIZE = 50
REQUEST_TIMEOUT_GET = 30
REQUEST_TIMEOUT_POST = 60
MAX_RETRIES = 3
BATCH_SLEEP_SECONDS = 0.5

DEFAULT_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
}


class EnsemblClient:
    def __init__(self) -> None:
        self.session = requests.Session()

    def _make_request(
        self,
        endpoint: str,
        method: str = "GET",
        json_data: Optional[Dict[str, Any]] = None,
    ) -> Optional[Any]:
        url = f"{BASE_URL}{endpoint}"
        for attempt in range(MAX_RETRIES):
            try:
                if method.upper() == "POST":
                    response = self.session.post(
                        url,
                        headers=DEFAULT_HEADERS,
                        json=json_data,
                        timeout=REQUEST_TIMEOUT_POST,
                    )
                else:
                    response = self.session.get(
                        url,
                        headers=DEFAULT_HEADERS,
                        timeout=REQUEST_TIMEOUT_GET,
                    )

                if response.status_code == 200:
                    return response.json()

                if response.status_code == 429 or 500 <= response.status_code < 600:
                    wait_time = int(response.headers.get("Retry-After", 2 ** attempt))
                    print(
                        f"Temporary Ensembl REST issue HTTP {response.status_code}. "
                        f"Waiting {wait_time}s...",
                        file=sys.stderr,
                    )
                    time.sleep(wait_time)
                    continue

                print(
                    f"HTTP error {response.status_code} for {url}: {response.text}",
                    file=sys.stderr,
                )
                return None

            except requests.exceptions.RequestException as exc:
                print(
                    f"Request error attempt {attempt + 1}/{MAX_RETRIES}: {exc}",
                    file=sys.stderr,
                )
                if attempt < MAX_RETRIES - 1:
                    time.sleep(2 ** attempt)
                else:
                    return None
        return None

    def get_gene_info(self, gene_id: str) -> Dict[str, Any]:
        data = self._make_request(f"/lookup/id/{gene_id}?expand=1")
        if not data:
            raise ValueError(
                f"Gene {gene_id} was not found, or Ensembl REST returned an error."
            )
        return data

    def get_sequences_batch(self, ids: List[str], seq_type: str) -> Dict[str, str]:
        unique_ids = sorted(set(ids))
        if not unique_ids:
            return {}

        all_results: Dict[str, str] = {}
        for i in range(0, len(unique_ids), BATCH_SIZE):
            batch = unique_ids[i : i + BATCH_SIZE]
            payload = {"ids": batch, "type": seq_type}

            print(
                f"  Downloading {seq_type} batch {i // BATCH_SIZE + 1} "
                f"({len(batch)} sequences)...",
                file=sys.stderr,
            )

            data = self._make_request("/sequence/id", method="POST", json_data=payload)
            if not data:
                continue
            if isinstance(data, dict):
                data = [data]

            for item in data:
                seq_id = item.get("id")
                sequence = item.get("seq", "")
                if seq_id:
                    all_results[seq_id] = sequence

            if i + BATCH_SIZE < len(unique_ids):
                time.sleep(BATCH_SLEEP_SECONDS)

        return all_results


