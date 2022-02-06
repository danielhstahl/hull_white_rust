window.BENCHMARK_DATA = {
  "lastUpdate": 1644171772771,
  "repoUrl": "https://github.com/danielhstahl/hull_white_rust",
  "entries": {
    "Benchmark": [
      {
        "commit": {
          "author": {
            "email": "danstahl1138@gmail.com",
            "name": "Daniel Stahl",
            "username": "danielhstahl"
          },
          "committer": {
            "email": "danstahl1138@gmail.com",
            "name": "Daniel Stahl",
            "username": "danielhstahl"
          },
          "distinct": true,
          "id": "14364a56cf21c4dc22bd06eed47a9f7d599355a9",
          "message": "updated actions",
          "timestamp": "2022-02-06T12:20:30-06:00",
          "tree_id": "25a17044edd1bb486fd1bf09811606b92395895f",
          "url": "https://github.com/danielhstahl/hull_white_rust/commit/14364a56cf21c4dc22bd06eed47a9f7d599355a9"
        },
        "date": 1644171771989,
        "tool": "cargo",
        "benches": [
          {
            "name": "bench_bond_now",
            "value": 31,
            "range": "± 4",
            "unit": "ns/iter"
          },
          {
            "name": "bench_bond_t",
            "value": 50,
            "range": "± 8",
            "unit": "ns/iter"
          },
          {
            "name": "bench_coupon_bond_now",
            "value": 146,
            "range": "± 19",
            "unit": "ns/iter"
          },
          {
            "name": "bench_coupon_bond_t",
            "value": 407,
            "range": "± 75",
            "unit": "ns/iter"
          },
          {
            "name": "bench_swap_rate",
            "value": 1595,
            "range": "± 251",
            "unit": "ns/iter"
          },
          {
            "name": "bench_swaption_american",
            "value": 35466367,
            "range": "± 5009158",
            "unit": "ns/iter"
          },
          {
            "name": "bench_swaption_european",
            "value": 22473,
            "range": "± 2355",
            "unit": "ns/iter"
          }
        ]
      }
    ]
  }
}