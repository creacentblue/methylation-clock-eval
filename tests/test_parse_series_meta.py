import gzip
import importlib

parse = importlib.import_module("03_parse")

FIXTURE = '''!Series_title\t"tiny synthetic fixture"
!Sample_title\t"age 67y 1001"\t"age 40y 1002"
!Sample_geo_accession\t"GSM1"\t"GSM2"
!Sample_source_name_ch1\t"whole blood"\t"whole blood"
!Sample_characteristics_ch1\t"age (y): 67"\t"age (y): 40"
!Sample_characteristics_ch1\t"gender: F"\t"gender: M"
!series_matrix_table_begin
"ID_REF"\t"GSM1"\t"GSM2"
"cg0001"\t0.1\t0.2
!series_matrix_table_end
'''


def test_series_meta_parses_characteristics(tmp_path):
    f = tmp_path / "tiny_series_matrix.txt.gz"
    with gzip.open(f, "wt") as fh:
        fh.write(FIXTURE)
    meta = parse.series_meta(f)
    assert list(meta.index) == ["GSM1", "GSM2"]
    assert meta.loc["GSM1", "age (y)"] == "67" and meta.loc["GSM2", "gender"] == "M"


def test_needed_set_empty_at_import():
    # importing must not require downloaded/processed data
    assert isinstance(parse.NEEDED, set)
