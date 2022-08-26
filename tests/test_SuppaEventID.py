from astk.event import SuppaEventID


def test_SE_SuppaEventID():
    pos_eid = "ENSMUSG00000025903.14;SE:chr1:4828649-4830268:4830315-4832311:+"
    neg_eid = "ENSMUSG00000025900.13;SE:chr1:4293012-4311270:4311433-4351910:-"

    pos_eidi = SuppaEventID(pos_eid)
    neg_eidi = SuppaEventID(neg_eid)

    assert pos_eidi.gene_id == "ENSMUSG00000025903.14"
    assert neg_eidi.gene_id == "ENSMUSG00000025900.13"
    assert pos_eidi.AS_type == "SE"
    assert neg_eidi.AS_type == "SE"
    assert pos_eidi.Chr == "chr1"
    assert neg_eidi.Chr == "chr1"
    assert pos_eidi.strand == "+"
    assert neg_eidi.strand == "-"

    test_pos_coor = (4828649, 4830268, 4830315, 4832311)
    test_neg_coor = (4293012, 4311270, 4311433, 4351910)
    assert pos_eidi.coordinates == test_pos_coor
    assert neg_eidi.coordinates == test_neg_coor
    assert pos_eidi.coordinate_str == "-".join(map(str, test_pos_coor))
    assert neg_eidi.coordinate_str == "-".join(map(str, test_neg_coor))
    assert pos_eidi.alter_element_coor == (4830268, 4830315)
    assert neg_eidi.alter_element_coor == (4311270, 4311433)
    assert pos_eidi.alter_element_len == 47
    assert neg_eidi.alter_element_len == 163


def test_A5_SuppaEventID():
    pos_eid = "ENSMUSG00000056763.16;A5:chr1:10039887-10045076:10039879-10045076:+"
    neg_eid = "ENSMUSG00000025917.9;A5:chr1:10035163-10037662:10035163-10037670:-"

    pos_eidi = SuppaEventID(pos_eid)
    neg_eidi = SuppaEventID(neg_eid)

    assert pos_eidi.gene_id == "ENSMUSG00000056763.16"
    assert neg_eidi.gene_id == "ENSMUSG00000025917.9"
    assert pos_eidi.AS_type == "A5"
    assert neg_eidi.AS_type == "A5"
    assert pos_eidi.Chr == "chr1"
    assert neg_eidi.Chr == "chr1"
    assert pos_eidi.strand == "+"
    assert neg_eidi.strand == "-"

    test_pos_coor = (10039879, 10039887, 10045076)
    test_neg_coor = (10035163, 10037662, 10037670)
    assert pos_eidi.coordinates == test_pos_coor
    assert neg_eidi.coordinates == test_neg_coor
    assert pos_eidi.coordinate_str == "10039887-10045076-10039879-10045076"
    assert neg_eidi.coordinate_str == "10035163-10037662-10035163-10037670"
    assert pos_eidi.alter_element_coor == (10039879, 10039887)
    assert neg_eidi.alter_element_coor == (10037662, 10037670)
    assert pos_eidi.alter_element_len == 8
    assert neg_eidi.alter_element_len == 8


def test_A3_SuppaEventID():
    pos_eid = "ENSMUSG00000005980.15;A3:chr16:4037174-4037308:4037174-4037357:+"
    neg_eid = "ENSMUSG00000098234.7;A3:chr1:9942155-9942323:9942149-9942323:-"
    pos_eidi = SuppaEventID(pos_eid)
    neg_eidi = SuppaEventID(neg_eid)


    assert pos_eidi.gene_id == "ENSMUSG00000005980.15"
    assert neg_eidi.gene_id == "ENSMUSG00000098234.7"
    assert pos_eidi.AS_type == "A3"
    assert neg_eidi.AS_type == "A3"
    assert pos_eidi.Chr == "chr16"
    assert neg_eidi.Chr == "chr1"
    assert pos_eidi.strand == "+"
    assert neg_eidi.strand == "-"

    test_pos_coor = (4037174, 4037308, 4037357)
    test_neg_coor = (9942149, 9942155, 9942323)
    assert pos_eidi.coordinates == test_pos_coor
    assert neg_eidi.coordinates == test_neg_coor
    assert pos_eidi.coordinate_str == "4037174-4037308-4037174-4037357"
    assert neg_eidi.coordinate_str == "9942155-9942323-9942149-9942323"
    assert pos_eidi.alter_element_coor == (4037308, 4037357)
    assert neg_eidi.alter_element_coor == (9942149, 9942155)
    assert pos_eidi.alter_element_len == 49
    assert neg_eidi.alter_element_len == 6


def test_AF_SuppaEventID():
    pos_eid = "ENSMUSG00000038822.15;AF:chr10:45577829:45578150-45588380:45578280:45578414-45588380:+"
    neg_eid = "ENSMUSG00000002459.17;AF:chr1:4923989-5019311:5019379:4923989-5070010:5070282:-"
    
    pos_eidi = SuppaEventID(pos_eid)
    neg_eidi = SuppaEventID(neg_eid)

    assert pos_eidi.gene_id == "ENSMUSG00000038822.15"
    assert neg_eidi.gene_id == "ENSMUSG00000002459.17"
    assert pos_eidi.AS_type == "AF"
    assert neg_eidi.AS_type == "AF"
    assert pos_eidi.Chr == "chr10"
    assert neg_eidi.Chr == "chr1"
    assert pos_eidi.strand == "+"
    assert neg_eidi.strand == "-"

    test_pos_coor = (45577829, 45578150, 45578280, 45578414, 45588380)
    test_neg_coor = (4923989, 5019311, 5019379, 5070010, 5070282)
    assert pos_eidi.coordinates == test_pos_coor
    assert neg_eidi.coordinates == test_neg_coor
    assert pos_eidi.coordinate_str == "45577829-45578150-45588380-45578280-45578414-45588380"
    assert neg_eidi.coordinate_str == "4923989-5019311-5019379-4923989-5070010-5070282"
    assert pos_eidi.alter_element_coor == (45577829, 45578150)
    assert neg_eidi.alter_element_coor == (5070010, 5070282)
    assert pos_eidi.alter_element_len == 321
    assert neg_eidi.alter_element_len == 272


def test_Al_SuppaEventID():
    pos_eid = "ENSMUSG00000030532.6;AL:chr7:80343956-80345177:80345315:80343956-80345648:80345765:+"
    neg_eid = "ENSMUSG00000025920.19;AL:chr1:16440442:16440458-16463035:16459715:16460419-16463035:-"

    pos_eidi = SuppaEventID(pos_eid)
    neg_eidi = SuppaEventID(neg_eid)

    assert pos_eidi.gene_id == "ENSMUSG00000030532.6"
    assert neg_eidi.gene_id == "ENSMUSG00000025920.19"
    assert pos_eidi.AS_type == "AL"
    assert neg_eidi.AS_type == "AL"
    assert pos_eidi.Chr == "chr7"
    assert neg_eidi.Chr == "chr1"
    assert pos_eidi.strand == "+"
    assert neg_eidi.strand == "-"

    test_pos_coor = (80343956, 80345177, 80345315, 80345648, 80345765)
    test_neg_coor = (16440442, 16440458, 16459715, 16460419, 16463035)
    assert pos_eidi.coordinates == test_pos_coor
    assert neg_eidi.coordinates == test_neg_coor
    assert pos_eidi.coordinate_str == "80343956-80345177-80345315-80343956-80345648-80345765"
    assert neg_eidi.coordinate_str == "16440442-16440458-16463035-16459715-16460419-16463035"
    assert pos_eidi.alter_element_coor == (80345648, 80345765)
    assert neg_eidi.alter_element_coor == (16440442, 16440458)
    assert pos_eidi.alter_element_len == 117
    assert neg_eidi.alter_element_len == 16


def test_MX_SuppaEventID():
    pos_eid = "ENSMUSG00000026020.9;MX:chr1:59704208-59704624:59704750-59706633:59704208-59705728:59705756-59706633:+"
    neg_eid = "ENSMUSG00000025968.16;MX:chr1:63170914-63172154:63172257-63176428:63170914-63175637:63175688-63176428:-"

    pos_eidi = SuppaEventID(pos_eid)
    neg_eidi = SuppaEventID(neg_eid)

    assert pos_eidi.gene_id == "ENSMUSG00000026020.9"
    assert neg_eidi.gene_id == "ENSMUSG00000025968.16"
    assert pos_eidi.AS_type == "MX"
    assert neg_eidi.AS_type == "MX"
    assert pos_eidi.Chr == "chr1"
    assert neg_eidi.Chr == "chr1"
    assert pos_eidi.strand == "+"
    assert neg_eidi.strand == "-"

    test_pos_coor = (59704208, 59704624, 59704750, 59705728, 59705756, 59706633)
    test_neg_coor = (63170914, 63172154, 63172257, 63175637, 63175688, 63176428)
    assert pos_eidi.coordinates == test_pos_coor
    assert neg_eidi.coordinates == test_neg_coor
    assert pos_eidi.coordinate_str == "59704208-59704624-59704750-59706633-59704208-59705728-59705756-59706633"
    assert neg_eidi.coordinate_str == "63170914-63172154-63172257-63176428-63170914-63175637-63175688-63176428"
    assert pos_eidi.alter_element_coor == (59704624, 59704750)
    assert neg_eidi.alter_element_coor == (63172154, 63172257)
    assert pos_eidi.alter_element_len == 126
    assert neg_eidi.alter_element_len == 103


def test_RI_SuppaEventID():
    pos_eid = "ENSMUSG00000026127.13;RI:chr1:34443769:34443923-34443997:34444091:+"
    neg_eid = "ENSMUSG00000026121.13;RI:chr1:36551726:36551838-36551927:36552149:-"

    pos_eidi = SuppaEventID(pos_eid)
    neg_eidi = SuppaEventID(neg_eid)

    assert pos_eidi.gene_id == "ENSMUSG00000026127.13"
    assert neg_eidi.gene_id == "ENSMUSG00000026121.13"
    assert pos_eidi.AS_type == "RI"
    assert neg_eidi.AS_type == "RI"
    assert pos_eidi.Chr == "chr1"
    assert neg_eidi.Chr == "chr1"
    assert pos_eidi.strand == "+"
    assert neg_eidi.strand == "-"

    test_pos_coor = (34443769, 34443923, 34443997, 34444091)
    test_neg_coor = (36551726, 36551838, 36551927, 36552149)
    assert pos_eidi.coordinates == test_pos_coor
    assert neg_eidi.coordinates == test_neg_coor
    assert pos_eidi.coordinate_str == "34443769-34443923-34443997-34444091"
    assert neg_eidi.coordinate_str == "36551726-36551838-36551927-36552149"
    assert pos_eidi.alter_element_coor == (34443923, 34443997)
    assert neg_eidi.alter_element_coor == (36551838, 36551927)
    assert pos_eidi.alter_element_len == 74
    assert neg_eidi.alter_element_len == 89
