import os
import unittest
from unittest.mock import patch

import pytest

from src import vcf_qc
from tests import TEST_DATA_DIR


class TestIntersectVcfWithBed(unittest.TestCase):
    def test_correct_unique_variants_output_when_intersected_against_bed(self):
        """
        Test that when we intersect against a bed file with overlapping
        regions we only retain unique variants (from the `-u` bedtools
        argument) and they are within the bed file regions
        """
        tmp_vcf = vcf_qc.intersect_vcf_with_bed(
            vcf=os.path.join(TEST_DATA_DIR, "test.vcf"),
            bed=os.path.join(TEST_DATA_DIR, "test.bed"),
        )

        intersected_variants = []

        with open(tmp_vcf) as fh:
            for line in fh.readlines():
                if not line.startswith("#"):
                    intersected_variants.append(line.split("\t")[:5])

        expected_variants = [
            ["1", "11119899", "rs7545802", "T", "C"],
            ["1", "17380497", "rs2746462", "G", "T"],
            ["1", "18807536", "rs7512414", "G", "T"],
        ]

        os.remove(tmp_vcf)

        with self.subTest("correct variants"):
            self.assertEqual(expected_variants, intersected_variants)

    @patch("src.vcf_qc.subprocess.run", wraps=vcf_qc.subprocess.run)
    def test_non_zero_exit_code_from_bedtools_raises_assertion_error(self, mock_run):

        mock_run.return_value.returncode = 1

        with pytest.raises(AssertionError):
            vcf_qc.intersect_vcf_with_bed(
                vcf=os.path.join(TEST_DATA_DIR, "test.vcf"),
                bed=os.path.join(TEST_DATA_DIR, "test.bed"),
            )


class TestIsAutosome(unittest.TestCase):
    def test_autosomes_correctly_return_true(self):
        autosomes = [
            "1",
            "chr1",
            "2",
            "chr2",
            "3",
            "chr3",
            "4",
            "chr4",
            "5",
            "chr5",
            "6",
            "chr6",
            "7",
            "chr7",
            "8",
            "chr8",
            "9",
            "chr9",
            "10",
            "chr10",
            "11",
            "chr11",
            "12",
            "chr12",
            "13",
            "chr13",
            "14",
            "chr14",
            "15",
            "chr15",
            "16",
            "chr16",
            "17",
            "chr17",
            "18",
            "chr18",
            "19",
            "chr19",
            "20",
            "chr20",
            "21",
            "chr21",
            "22",
            "chr22",
        ]

        for chrom in autosomes:
            with self.subTest("automsomes return True"):
                assert vcf_qc.is_autosome(chrom)

    def test_non_autosomes_correctly_return_false(self):
        non_autosomes = ["X", "chrX", "Y", "chrY", "GL000192.1"]

        for chrom in non_autosomes:
            with self.subTest("automsomes return True"):
                assert not vcf_qc.is_autosome(chrom)


class TestGetHetHomCounts(unittest.TestCase):
    """
    Tests for calculating the AAF of all het and hom variants in the vcf.

    We are calculating AFF per variant as the AD of all the non reference
    alleles / the sum of all ADs at that position. We are not using the
    DP to ensure non informative reads do not contribute to the depth.
    """

    def test_aaf_calculated_correctly_for_hets_and_homs(self):
        """
        Test using our minimal `test.vcf` with 2 hets and 2 homs that
        the returned values are correct.

        Expected variant values:
            GT      ADs         total AD    expected AAF
            (1, 1) (0, 83)      83          1.0
            (1, 1) (0, 1367)    1367        1.0
            (0, 1) (158, 155)   313         0.4952076677316294
            (0, 1) (60, 50)     110         0.45454545454545453

        """
        calculated_values = vcf_qc.get_het_hom_counts(
            os.path.join(TEST_DATA_DIR, "test.vcf")
        )

        expected_values = {
            "het": [0.45454545454545453, 0.4952076677316294],
            "hom": [1.0, 1.0],
            "x_het": [],
            "x_hom": [],
        }

        self.assertEqual(expected_values, calculated_values)

    def test_x_chrom_variants_additionally_returned_separate(self):
        """
        Test that X chromosome variants are correctly returned in a
        separate keys but also contribute towards normal het hom
        counts

        Expected chromosome 1 variant values:
            GT      ADs         total AD    expected AAF
            (1, 1) (0, 83)      83          1.0
            (1, 1) (0, 1367)    1367        1.0
            (0, 1) (158, 155)   313         0.4952076677316294
            (0, 1) (60, 50)     110         0.45454545454545453

        Expected chromsome X variants
            GT      ADs         total AD    expected AAF
            (0, 1)  (1, 4)      5           0.8
            (0, 1)  (2, 3)      5           0.6
            (1, 1)  (0, 2)      2           1.0
            (1, 1)  (0, 4)      4           1.0
        """
        calculated_values = vcf_qc.get_het_hom_counts(
            os.path.join(TEST_DATA_DIR, "test_w_x.vcf")
        )

        expected_values = {
            "het": [0.45454545454545453, 0.4952076677316294],
            "hom": [1.0, 1.0],
            "x_het": [0.6, 0.8],
            "x_hom": [1.0, 1.0],
        }

        self.assertEqual(calculated_values, expected_values)

    def test_reference_calls_are_skipped(self):
        """
        Test that if any reference calls are present that these are
        skipped and not included in the calculations as this suggests
        a gVCF has been passed.

        The `test_w_ref.vcf` has only reference variants, therefore we
        should get back empty lists
        """
        calculated_values = vcf_qc.get_het_hom_counts(
            os.path.join(TEST_DATA_DIR, "test_w_ref.vcf")
        )

        expected_values = {"het": [], "hom": [], "x_het": [], "x_hom": []}

        self.assertEqual(calculated_values, expected_values)

    def test_multi_sample_vcf_raises_assertion_error(self):

        match = "more than one sample present in vcf"
        with pytest.raises(AssertionError, match=match):
            calculated_values = vcf_qc.get_het_hom_counts(
                os.path.join(TEST_DATA_DIR, "multi_sample.vcf")
            )


class TestCalculateRatios(unittest.TestCase):

    def test_het_and_hom_ratios_calculated_correctly(self):
        """
        Test that the mean het, mean hom and het:hom ratios are correct.

        Mean het and hom are just the mean of AAFs from get_het_hom_counts,
        and het:hom ratio is the ratio of total het to total non-ref hom
        variants
        """
        counts = {
            "het": [0.5, 0.45, 0.55],
            "hom": [1.0, 1.0, 1.0, 1.0, 0.95],
            "x_het": [],
            "x_hom": [],
        }

        expected_ratios = {
            "mean_het": "0.5000",
            "mean_hom": "0.9900",
            "het_hom_ratio": "0.6000",
            "x_het_hom_ratio": None,
        }

        calculated_ratios = vcf_qc.calculate_ratios(counts)

        self.assertEqual(expected_ratios, calculated_ratios)

    def test_x_het_hom_calculated_correctly(self):
        counts = {
            "het": [0.5, 0.45, 0.55],
            "hom": [1.0, 1.0, 1.0, 1.0, 0.95],
            "x_het": [0.45, 0.55],
            "x_hom": [1.0, 0.98],
        }

        calculated_ratios = vcf_qc.calculate_ratios(counts)

        self.assertEqual(calculated_ratios["x_het_hom_ratio"], "1.0000")

    def test_none_values_returned_if_missing_het_or_hom_counts(self):
        empty_ratios = {
            "mean_het": None,
            "mean_hom": None,
            "het_hom_ratio": None,
            "x_het_hom_ratio": None,
        }

        with self.subTest("missing het"):
            self.assertEqual(
                vcf_qc.calculate_ratios({"het": [], "hom": ["0.5", "0.5"]}),
                empty_ratios,
            )

        with self.subTest("missing hom"):
            self.assertEqual(
                vcf_qc.calculate_ratios({"het": ["1.0", "1.0"], "hom": []}),
                empty_ratios,
            )

        with self.subTest("missing het and hom"):
            self.assertEqual(
                vcf_qc.calculate_ratios({"het": [], "hom": []}), empty_ratios
            )


@patch("src.vcf_qc.dxpy.describe", return_value={"name": "sample1.vcf"})
@patch("src.vcf_qc.dxpy.bindings.dxfile_functions.download_dxfile")
class TestDownloadInputFile(unittest.TestCase):

    def test_remote_filename_returned(self, mock_download, mock_describe):
        self.assertEqual(vcf_qc.download_input_file("file-xxx"), "sample1.vcf")

    def test_download_called_with_params(self, mock_download, mock_describe):
        vcf_qc.download_input_file("file-xxx")

        self.assertEqual(
            mock_download.call_args[1], {"dxid": "file-xxx", "filename": "sample1.vcf"}
        )


class TestWriteOutputFile(unittest.TestCase):
    def test_output_file_contents_correctly_written(self):
        ratios = {
            "mean_het": "0.5000",
            "mean_hom": "0.9900",
            "het_hom_ratio": "0.6000",
            "x_het_hom_ratio": None,
        }

        vcf_qc.write_output_file(outfile="test.vcf.qc", ratios=ratios)

        with open("test.vcf.qc", "r") as fh:
            written_contents = fh.readlines()

        expected_contents = [
            "mean_het\tmean_hom\thet_hom_ratio\tx_het_hom_ratio\n",
            "0.5000\t0.9900\t0.6000\tNone\n",
        ]

        os.remove("test.vcf.qc")

        self.assertEqual(written_contents, expected_contents)


@patch("src.vcf_qc.dxpy.upload_local_file")
@patch("src.vcf_qc.dxpy.bindings.dxjob.DXJob")
class TestUploadOutputFile(unittest.TestCase):
    def test_upload_params_correct(self, mock_job, mock_upload):
        mock_job.return_value.describe.return_value = {"folder": "/output/sub_folder"}

        vcf_qc.upload_output_file("test.vcf.qc")

        expected_args = {
            "filename": "test.vcf.qc",
            "folder": "/output/sub_folder",
            "wait_on_close": True,
        }

        self.assertEqual(mock_upload.call_args[1], expected_args)

    def test_expected_dxlink_format_returned(self, mock_job, mock_upload):
        mock_upload.return_value = "file-xxx"

        self.assertEqual(
            vcf_qc.upload_output_file("test.vcf.qc"),
            {"output_file": {"$dnanexus_link": "file-xxx"}},
        )
