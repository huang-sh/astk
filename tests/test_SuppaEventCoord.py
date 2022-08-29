from astk.utils import SuppaEventCoord


def test_SE_SuppaEventCoord():
    pos_eid = "ENSMUSG00000025903.14;SE:chr1:4828649-4830268:4830315-4832311:+"
    neg_eid = "ENSMUSG00000025900.13;SE:chr1:4293012-4311270:4311433-4351910:-"

    pos_eidi = SuppaEventCoord(pos_eid)
    neg_eidi = SuppaEventCoord(neg_eid)

    assert pos_eidi.get_alter(False)  == (4830268, 4830315)
    assert pos_eidi.get_alter()  == (4830267, 4830315)

    assert neg_eidi.get_alter(False) == (4311270, 4311433)
    assert neg_eidi.get_alter() == (4311269, 4311433)


    assert pos_eidi.get_ss_range(1, 4)  == (4828648, 4832311)
    assert pos_eidi.get_ss_range(2, 3, based0=False) == (4830268, 4830315)
    assert pos_eidi.get_ss_range(1, 3, True, False) == (4828649, 4830315)
    pos_kw_res = pos_eidi.get_ss_range(1, 2, based0=False, s_ups=50, s_dws=10,e_ups=10)
    assert pos_kw_res == (4828649-50+10, 4830268-10)

    assert neg_eidi.get_ss_range(2, 4)  == (4311269, 4351910)
    assert neg_eidi.get_ss_range(2, 3, based0=False) == (4311270, 4311433)
    assert neg_eidi.get_ss_range(1, 3, True, False) == (4311270, 4351910)
    neg_kw_res = neg_eidi.get_ss_range(1, 2, strand_sp=True, based0=False, s_ups=50, s_dws=10,e_ups=10)
    assert neg_kw_res == (4311433+10, 4351910+50-10)

    assert pos_eidi.get_splicescore_flank(1)  == (4828649-2-1, 4828649+6)
    assert pos_eidi.get_splicescore_flank(1, based0=False) == (4828649-2, 4828649+6)
    assert pos_eidi.get_splicescore_flank(2) == (4830268-20-1, 4830268+2)
    assert pos_eidi.get_splicescore_flank(2, False) == (4830268-20, 4830268+2)
    assert pos_eidi.get_splicescore_flank(3) == (4830315-2-1, 4830315+6)
    assert pos_eidi.get_splicescore_flank(3, False) == (4830315-2, 4830315+6)
    assert pos_eidi.get_splicescore_flank(4) == (4832311-20-1, 4832311+2)
    assert pos_eidi.get_splicescore_flank(4, False) == (4832311-20, 4832311+2) 

    assert neg_eidi.get_splicescore_flank(1)  == (4351910-6-1, 4351910+2)
    assert neg_eidi.get_splicescore_flank(1, based0=False) == (4351910-6, 4351910+2)
    assert neg_eidi.get_splicescore_flank(2) == (4311433-2-1, 4311433+20)
    assert neg_eidi.get_splicescore_flank(2, False) == (4311433-2, 4311433+20)
    assert neg_eidi.get_splicescore_flank(3) == (4311270-6-1, 4311270+2)
    assert neg_eidi.get_splicescore_flank(3, False) == (4311270-6, 4311270+2)
    assert neg_eidi.get_splicescore_flank(4) == (4293012-2-1, 4293012+20)
    assert neg_eidi.get_splicescore_flank(4, False) == (4293012-2, 4293012+20)


    assert pos_eidi.get_ss_flank(1,sss=True)  == (4828649-2-1, 4828649+6)
    assert pos_eidi.get_ss_flank(2,sss=True) == (4830268-20-1, 4830268+2)
    assert pos_eidi.get_ss_flank(3,sss=True) == (4830315-2-1, 4830315+6)
    assert pos_eidi.get_ss_flank(4,sss=True) == (4832311-20-1, 4832311+2)
    assert neg_eidi.get_ss_flank(1,sss=True)  == (4351910-6-1, 4351910+2)
    assert neg_eidi.get_ss_flank(2,sss=True) == (4311433-2-1, 4311433+20)
    assert neg_eidi.get_ss_flank(3,sss=True) == (4311270-6-1, 4311270+2)
    assert neg_eidi.get_ss_flank(4,sss=True) == (4293012-2-1, 4293012+20)

    assert pos_eidi.get_ss_flank_bed(1)  == ("chr1", 4828649-150-1, 4828649+150, pos_eid, 0, "+")
    assert pos_eidi.get_ss_flank_bed(1, sss=True)  == ("chr1",4828649-2-1, 4828649+6, pos_eid, 0, "+")
    assert neg_eidi.get_ss_flank_bed(1,sss=True)  == ("chr1" , 4351910-6-1, 4351910+2, neg_eid, 0, "-")

    assert pos_eidi.get_ss_flank(1,)  == (4828649-150-1, 4828649+150)
    assert pos_eidi.get_ss_flank(2, ups_w=100) == (4830268-100-1, 4830268+150)
    assert pos_eidi.get_ss_flank(3,dws_w=100,based0=False) == (4830315-150, 4830315+100)
    assert pos_eidi.get_ss_flank(4,ups_w=100, dws_w=100,based0=False) == (4832311-100, 4832311+100)

    assert neg_eidi.get_ss_flank(1,)  == (4293012-150-1, 4293012+150)
    assert neg_eidi.get_ss_flank(2, ups_w=100) == (4311270-100-1, 4311270+150)
    assert neg_eidi.get_ss_flank(3,dws_w=100,based0=False) == (4311433-150, 4311433+100)
    assert neg_eidi.get_ss_flank(4,ups_w=100, dws_w=100,based0=False) == (4351910-100, 4351910+100)

    assert neg_eidi.get_ss_flank(1,strand_sp=True)  == (4351910-150-1, 4351910+150)
    assert neg_eidi.get_ss_flank(2, ups_w=100,strand_sp=True) == (4311433-150-1, 4311433+100)
    assert neg_eidi.get_ss_flank(3,dws_w=100,strand_sp=True,based0=False) == (4311270-100, 4311270+150)
    assert neg_eidi.get_ss_flank(4,ups_w=100, dws_w=110,strand_sp=True,based0=False) == (4293012-110, 4293012+100)


def test_A5_SuppaEventCoord():
    pos_eid = "ENSMUSG00000056763.16;A5:chr1:10039887-10045076:10039879-10045076:+"
    neg_eid = "ENSMUSG00000025917.9;A5:chr1:10035163-10037662:10035163-10037670:-"

    pos_eidi = SuppaEventCoord(pos_eid)
    neg_eidi = SuppaEventCoord(neg_eid)

    assert pos_eidi.get_alter(False)  == (10039879, 10039887)
    assert pos_eidi.get_alter()  == (10039879-1, 10039887)
    assert neg_eidi.get_alter(False) == (10037662, 10037670)
    assert neg_eidi.get_alter() ==  (10037662-1, 10037670)


    assert pos_eidi.get_ss_range(1, 2)  == (10039879-1, 10039887)
    assert pos_eidi.get_ss_range(2, 3, based0=False) == (10039887, 10045076)
    assert pos_eidi.get_ss_range(2, 3, True, False) == (10039887, 10045076)

    assert neg_eidi.get_ss_range(1, 2)  == (10035163-1, 10037662)
    assert neg_eidi.get_ss_range(2, 3, based0=False) == (10037662, 10037670)
    assert neg_eidi.get_ss_range(2, 3, True, False) == (10035163, 10037662)

    assert pos_eidi.get_splicescore_flank(1)  == (10039879-2-1, 10039879+6)
    assert pos_eidi.get_splicescore_flank(1, based0=False) == (10039879-2, 10039879+6)
    assert pos_eidi.get_splicescore_flank(2) == (10039887-2-1, 10039887+6)
    assert pos_eidi.get_splicescore_flank(2, False) == (10039887-2, 10039887+6)
    assert pos_eidi.get_splicescore_flank(3) == (10045076-20-1, 10045076+2)
    assert pos_eidi.get_splicescore_flank(3, False) == (10045076-20, 10045076+2)


    assert neg_eidi.get_splicescore_flank(1)  == (10037670-6-1, 10037670+2)
    assert neg_eidi.get_splicescore_flank(1, based0=False) == (10037670-6, 10037670+2)
    assert neg_eidi.get_splicescore_flank(2) == (10037662-6-1, 10037662+2)
    assert neg_eidi.get_splicescore_flank(2, False) == (10037662-6, 10037662+2)
    assert neg_eidi.get_splicescore_flank(3) == (10035163-2-1, 10035163+20)
    assert neg_eidi.get_splicescore_flank(3, False) == (10035163-2, 10035163+20)


    assert pos_eidi.get_ss_flank(1,sss=True)  == (10039879-2-1, 10039879+6)
    assert pos_eidi.get_ss_flank(2,sss=True) == (10039887-2-1, 10039887+6)
    assert pos_eidi.get_ss_flank(3,sss=True) == (10045076-20-1, 10045076+2)
    assert neg_eidi.get_ss_flank(1,sss=True)  == (10037670-6-1, 10037670+2)
    assert neg_eidi.get_ss_flank(2,sss=True) == (10037662-6-1, 10037662+2)
    assert neg_eidi.get_ss_flank(3,sss=True) == (10035163-2-1, 10035163+20)


    assert pos_eidi.get_ss_flank(1,)  == (10039879-150-1, 10039879+150)
    assert pos_eidi.get_ss_flank(2, ups_w=100) == (10039887-100-1, 10039887+150)
    assert pos_eidi.get_ss_flank(3,dws_w=100,based0=False) == (10045076-150, 10045076+100)

    assert neg_eidi.get_ss_flank(1,)  == (10035163-150-1, 10035163+150)
    assert neg_eidi.get_ss_flank(2, ups_w=100) == (10037662-100-1, 10037662+150)
    assert neg_eidi.get_ss_flank(3,dws_w=100,based0=False) == (10037670-150, 10037670+100)

    assert neg_eidi.get_ss_flank(1,strand_sp=True)  == (10037670-150-1, 10037670+150)
    assert neg_eidi.get_ss_flank(2, ups_w=100,strand_sp=True) == (10037662-150-1, 10037662+100)
    assert neg_eidi.get_ss_flank(3,dws_w=100,strand_sp=True,based0=False) == (10035163-100, 10035163+150)
 

def test_A3_SuppaEventCoord():
    pos_eid = "ENSMUSG00000005980.15;A3:chr16:4037174-4037308:4037174-4037357:+"
    neg_eid = "ENSMUSG00000098234.7;A3:chr1:9942155-9942323:9942149-9942323:-"
    pos_eidi = SuppaEventCoord(pos_eid)
    neg_eidi = SuppaEventCoord(neg_eid)

    assert pos_eidi.get_alter(False)  == (4037308, 4037357)
    assert pos_eidi.get_alter()  == (4037308-1, 4037357)
    assert neg_eidi.get_alter(False) == (9942149, 9942155)
    assert neg_eidi.get_alter() ==  (9942149-1, 9942155)

    assert pos_eidi.get_ss_range(1, 2)  == (4037174-1, 4037308)
    assert pos_eidi.get_ss_range(2, 3, based0=False) == (4037308, 4037357)
    assert pos_eidi.get_ss_range(2, 3, True, False) == (4037308, 4037357)

    assert neg_eidi.get_ss_range(1, 2)  == (9942149-1, 9942155)
    assert neg_eidi.get_ss_range(2, 3, based0=False) == (9942155, 9942323)
    assert neg_eidi.get_ss_range(2, 3, True, False) == (9942149, 9942155)

    assert pos_eidi.get_splicescore_flank(1)  == (4037174-2-1, 4037174+6)
    assert pos_eidi.get_splicescore_flank(1, based0=False) == (4037174-2, 4037174+6)
    assert pos_eidi.get_splicescore_flank(2) == (4037308-20-1, 4037308+2)
    assert pos_eidi.get_splicescore_flank(2, False) == (4037308-20, 4037308+2)
    assert pos_eidi.get_splicescore_flank(3) == (4037357-20-1, 4037357+2)
    assert pos_eidi.get_splicescore_flank(3, False) == (4037357-20, 4037357+2)

    assert neg_eidi.get_splicescore_flank(1)  == (9942323-6-1, 9942323+2)
    assert neg_eidi.get_splicescore_flank(2) == (9942155-2-1, 9942155+20)
    assert neg_eidi.get_splicescore_flank(3) == (9942149-2-1, 9942149+20)

    assert pos_eidi.get_ss_flank(1,sss=True) == (4037174-2-1, 4037174+6)
    assert neg_eidi.get_ss_flank(1,sss=True) == (9942323-6-1, 9942323+2)


    assert pos_eidi.get_ss_flank(1,)  == (4037174-150-1, 4037174+150)
    assert pos_eidi.get_ss_flank(2, ups_w=100) == (4037308-100-1, 4037308+150)
    assert pos_eidi.get_ss_flank(3,dws_w=100,based0=False) == (4037357-150, 4037357+100)

    assert neg_eidi.get_ss_flank(1,)  == (9942149-150-1, 9942149+150)
    assert neg_eidi.get_ss_flank(2, ups_w=100) == (9942155-100-1, 9942155+150)
    assert neg_eidi.get_ss_flank(3,dws_w=100,based0=False) == (9942323-150, 9942323+100)

    assert neg_eidi.get_ss_flank(1,strand_sp=True)  == (9942323-150-1, 9942323+150)
    assert neg_eidi.get_ss_flank(2, ups_w=100,strand_sp=True) == (9942155-150-1, 9942155+100)
    assert neg_eidi.get_ss_flank(3,dws_w=100,strand_sp=True,based0=False) == (9942149-100, 9942149+150)


def test_AF_SuppaEventCoord():
    pos_eid = "ENSMUSG00000038822.15;AF:chr10:45577829:45578150-45588380:45578280:45578414-45588380:+"
    neg_eid = "ENSMUSG00000002459.17;AF:chr1:4923989-5019311:5019379:4923989-5070010:5070282:-"
    
    pos_eidi = SuppaEventCoord(pos_eid)
    neg_eidi = SuppaEventCoord(neg_eid)


    assert pos_eidi.get_alter(False)  == (45577829, 45578150)
    assert pos_eidi.get_alter()  == (45577829-1, 45578150)

    assert neg_eidi.get_alter(False) == (5070010, 5070282)
    assert neg_eidi.get_alter() == (5070010-1, 5070282)

    assert pos_eidi.get_ss_range(1, 4)  == (45577829-1, 45578414)
    assert neg_eidi.get_ss_range(2, 4)  == (5019311-1, 5070010)
    assert neg_eidi.get_ss_range(2, 5, True, False) == (4923989, 5070010)

    assert pos_eidi.get_splicescore_flank(1)  == (45577828, 45577829)
    assert pos_eidi.get_splicescore_flank(1, False) == (45577829, 45577829)
    assert pos_eidi.get_splicescore_flank(2) == (45578150-2-1, 45578150+6)
    assert pos_eidi.get_splicescore_flank(2, False) == (45578150-2, 45578150+6)
    assert pos_eidi.get_splicescore_flank(3) == (45578280-20-1, 45578280+2)
    assert pos_eidi.get_splicescore_flank(3, False) == (45578280-20, 45578280+2)
    assert pos_eidi.get_splicescore_flank(4) == (45578414-2-1, 45578414+6)
    assert pos_eidi.get_splicescore_flank(4, False) == (45578414-2, 45578414+6)
    assert pos_eidi.get_splicescore_flank(5) == (45588380-20-1, 45588380+2)
    assert pos_eidi.get_splicescore_flank(5, False) == (45588380-20, 45588380+2)

    assert neg_eidi.get_splicescore_flank(1)  == (5070282-1, 5070282)
    assert neg_eidi.get_splicescore_flank(1, False) == (5070282, 5070282)
    assert neg_eidi.get_splicescore_flank(2) == (5070010-6-1, 5070010+2)
    assert neg_eidi.get_splicescore_flank(2, False) == (5070010-6, 5070010+2)
    assert neg_eidi.get_splicescore_flank(3) == (5019379-2-1, 5019379+20)
    assert neg_eidi.get_splicescore_flank(3, False) == (5019379-2, 5019379+20)
    assert neg_eidi.get_splicescore_flank(4) == (5019311-6-1, 5019311+2)
    assert neg_eidi.get_splicescore_flank(4, False) == (5019311-6, 5019311+2)
    assert neg_eidi.get_splicescore_flank(5) == (4923989-2-1, 4923989+20)
    assert neg_eidi.get_splicescore_flank(5, False) == (4923989-2, 4923989+20)


    assert pos_eidi.get_ss_flank(1,sss=True)  == (45577829-1, 45577829)
    assert pos_eidi.get_ss_flank(2,sss=True) == (45578150-2-1, 45578150+6)
    assert pos_eidi.get_ss_flank(3,sss=True) == (45578280-20-1, 45578280+2)
    assert pos_eidi.get_ss_flank(4,sss=True) == (45578414-2-1, 45578414+6)
    assert neg_eidi.get_ss_flank(1,sss=True)  == (5070282-1, 5070282)
    assert neg_eidi.get_ss_flank(2,sss=True) == (5070010-6-1, 5070010+2)
    assert neg_eidi.get_ss_flank(3,sss=True) == (5019379-2-1, 5019379+20)
    assert neg_eidi.get_ss_flank(4,sss=True) == (5019311-6-1, 5019311+2)

    assert pos_eidi.get_ss_flank(1,)  == (45577829-150-1, 45577829+150)
    assert pos_eidi.get_ss_flank(2, ups_w=100) == (45578150-100-1, 45578150+150)
    assert pos_eidi.get_ss_flank(3,dws_w=100,based0=False) == (45578280-150, 45578280+100)
    assert pos_eidi.get_ss_flank(4,ups_w=100, dws_w=100,based0=False) == (45578414-100, 45578414+100)
    assert pos_eidi.get_ss_flank(5,ups_w=120, dws_w=100,based0=False) == (45588380-120, 45588380+100)

    assert neg_eidi.get_ss_flank(1,)  == (4923989-150-1, 4923989+150)
    assert neg_eidi.get_ss_flank(2, ups_w=100) == (5019311-100-1, 5019311+150)
    assert neg_eidi.get_ss_flank(3,dws_w=100,based0=False) == (5019379-150, 5019379+100)
    assert neg_eidi.get_ss_flank(4,ups_w=100, dws_w=100,based0=False) == (5070010-100, 5070010+100)

    assert neg_eidi.get_ss_flank(1,strand_sp=True)  == (5070282-150-1, 5070282+150)
    assert neg_eidi.get_ss_flank(2, ups_w=100,strand_sp=True) == (5070010-150-1, 5070010+100)
    assert neg_eidi.get_ss_flank(3,dws_w=100,strand_sp=True,based0=False) == (5019379-100, 5019379+150)
    assert neg_eidi.get_ss_flank(4,ups_w=100, dws_w=110,strand_sp=True,based0=False) == (5019311-110, 5019311+100)


def test_AL_SuppaEventCoord():
    pos_eid = "ENSMUSG00000030532.6;AL:chr7:80343956-80345177:80345315:80343956-80345648:80345765:+"
    neg_eid = "ENSMUSG00000025920.19;AL:chr1:16440442:16440458-16463035:16459715:16460419-16463035:-"

    pos_eidi = SuppaEventCoord(pos_eid)
    neg_eidi = SuppaEventCoord(neg_eid)

    assert pos_eidi.get_alter(False)  == (80345648, 80345765)
    assert pos_eidi.get_alter()  == (80345648-1, 80345765)

    assert neg_eidi.get_alter(False) == (16440442, 16440458)
    assert neg_eidi.get_alter() == (16440442-1, 16440458)

    assert pos_eidi.get_ss_range(1, 4)  == (80343956-1, 80345648)
    assert neg_eidi.get_ss_range(2, 4)  == (16440458-1, 16460419)
    assert neg_eidi.get_ss_range(2, 5, True, False) == (16440442, 16460419)

    assert pos_eidi.get_splicescore_flank(1)  == (80343956-2-1, 80343956+6)
    assert pos_eidi.get_splicescore_flank(1, False) == (80343956-2, 80343956+6)
    assert pos_eidi.get_splicescore_flank(2) == (80345177-20-1, 80345177+2)
    assert pos_eidi.get_splicescore_flank(2, False) == (80345177-20, 80345177+2)
    assert pos_eidi.get_splicescore_flank(3) == (80345315-2-1, 80345315+6)
    assert pos_eidi.get_splicescore_flank(3, False) == (80345315-2, 80345315+6)
    assert pos_eidi.get_splicescore_flank(4) == (80345648-20-1, 80345648+2)
    assert pos_eidi.get_splicescore_flank(4, False) == (80345648-20, 80345648+2)
    assert pos_eidi.get_splicescore_flank(5) == (80345765-1, 80345765)
    assert pos_eidi.get_splicescore_flank(5, False) == (80345765, 80345765)

    assert neg_eidi.get_splicescore_flank(1)  == (16463035-6-1, 16463035+2)
    assert neg_eidi.get_splicescore_flank(1, False) == (16463035-6, 16463035+2)
    assert neg_eidi.get_splicescore_flank(2) == (16460419-2-1, 16460419+20)
    assert neg_eidi.get_splicescore_flank(2, False) == (16460419-2, 16460419+20)
    assert neg_eidi.get_splicescore_flank(3) == (16459715-6-1, 16459715+2)
    assert neg_eidi.get_splicescore_flank(3, False) == (16459715-6, 16459715+2)
    assert neg_eidi.get_splicescore_flank(4) == (16440458-2-1, 16440458+20)
    assert neg_eidi.get_splicescore_flank(4, False) == (16440458-2, 16440458+20)
    assert neg_eidi.get_splicescore_flank(5) == (16440442-1, 16440442)
    assert neg_eidi.get_splicescore_flank(5, False) == (16440442, 16440442)


    assert pos_eidi.get_ss_flank(1,sss=True)  == (80343956-2-1, 80343956+6)
    assert pos_eidi.get_ss_flank(2,sss=True) == (80345177-20-1, 80345177+2)
    assert pos_eidi.get_ss_flank(3,sss=True) == (80345315-2-1, 80345315+6)
    assert pos_eidi.get_ss_flank(4,sss=True) ==(80345648-20-1, 80345648+2)
    assert neg_eidi.get_ss_flank(1,sss=True)  == (16463035-6-1, 16463035+2)
    assert neg_eidi.get_ss_flank(2,sss=True) == (16460419-2-1, 16460419+20)
    assert neg_eidi.get_ss_flank(3,sss=True) == (16459715-6-1, 16459715+2)
    assert neg_eidi.get_ss_flank(4,sss=True) == (16440458-2-1, 16440458+20)

    assert pos_eidi.get_ss_flank(1,)  == (80343956-150-1, 80343956+150)
    assert pos_eidi.get_ss_flank(2, ups_w=100) == (80345177-100-1, 80345177+150)
    assert pos_eidi.get_ss_flank(3,dws_w=100,based0=False) == (80345315-150, 80345315+100)
    assert pos_eidi.get_ss_flank(4,ups_w=100, dws_w=100,based0=False) == (80345648-100, 80345648+100)
    assert pos_eidi.get_ss_flank(5,ups_w=120, dws_w=100,based0=False) == (80345765-120, 80345765+100)

    assert neg_eidi.get_ss_flank(1,)  == (16440442-150-1, 16440442+150)
    assert neg_eidi.get_ss_flank(2, ups_w=100) == (16440458-100-1, 16440458+150)
    assert neg_eidi.get_ss_flank(3,dws_w=100,based0=False) == (16459715-150, 16459715+100)
    assert neg_eidi.get_ss_flank(4,ups_w=100, dws_w=100,based0=False) == (16460419-100, 16460419+100)

    assert neg_eidi.get_ss_flank(1,strand_sp=True)  == (16463035-150-1, 16463035+150)
    assert neg_eidi.get_ss_flank(2, ups_w=100,strand_sp=True) == (16460419-150-1, 16460419+100)
    assert neg_eidi.get_ss_flank(3,dws_w=100,strand_sp=True,based0=False) == (16459715-100, 16459715+150)
    assert neg_eidi.get_ss_flank(4,ups_w=100, dws_w=110,strand_sp=True,based0=False) == (16440458-110, 16440458+100)


def test_MX_SuppaEventCoord():
    pos_eid = "ENSMUSG00000026020.9;MX:chr1:59704208-59704624:59704750-59706633:59704208-59705728:59705756-59706633:+"
    neg_eid = "ENSMUSG00000025968.16;MX:chr1:63170914-63172154:63172257-63176428:63170914-63175637:63175688-63176428:-"

    pos_eidi = SuppaEventCoord(pos_eid)
    neg_eidi = SuppaEventCoord(neg_eid)


    assert pos_eidi.get_ss_range(1, 6)  == (59704208-1, 59706633)
    assert neg_eidi.get_ss_range(2, 4)  == (63172154-1, 63175637)
    assert neg_eidi.get_ss_range(2, 5, True, False) == (63172154, 63175688)

    assert pos_eidi.get_splicescore_flank(1)  == (59704208-2-1, 59704208+6)
    assert pos_eidi.get_splicescore_flank(1, False) == (59704208-2, 59704208+6)
    assert pos_eidi.get_splicescore_flank(2) == (59704624-20-1, 59704624+2)
    assert pos_eidi.get_splicescore_flank(2, False) == (59704624-20, 59704624+2)
    assert pos_eidi.get_splicescore_flank(3) == (59704750-2-1, 59704750+6)
    assert pos_eidi.get_splicescore_flank(3, False) == (59704750-2, 59704750+6)
    assert pos_eidi.get_splicescore_flank(4) == (59705728-20-1, 59705728+2)
    assert pos_eidi.get_splicescore_flank(4, False) == (59705728-20, 59705728+2)
    assert pos_eidi.get_splicescore_flank(5) == (59705756-2-1, 59705756+6)
    assert pos_eidi.get_splicescore_flank(5, False) == (59705756-2, 59705756+6)
    assert pos_eidi.get_splicescore_flank(6) == (59706633-20-1, 59706633+2)
    assert pos_eidi.get_splicescore_flank(6, False) == (59706633-20, 59706633+2)

    assert neg_eidi.get_splicescore_flank(1)  == (63176428-6-1, 63176428+2)
    assert neg_eidi.get_splicescore_flank(1, False) == (63176428-6, 63176428+2)
    assert neg_eidi.get_splicescore_flank(2) == (63175688-2-1, 63175688+20)
    assert neg_eidi.get_splicescore_flank(2, False) == (63175688-2, 63175688+20)
    assert neg_eidi.get_splicescore_flank(3) == (63175637-6-1, 63175637+2)
    assert neg_eidi.get_splicescore_flank(3, False) == (63175637-6, 63175637+2)
    assert neg_eidi.get_splicescore_flank(4) == (63172257-2-1, 63172257+20)
    assert neg_eidi.get_splicescore_flank(4, False) == (63172257-2, 63172257+20)
    assert neg_eidi.get_splicescore_flank(5) == (63172154-6-1, 63172154+2)
    assert neg_eidi.get_splicescore_flank(5, False) == (63172154-6, 63172154+2)
    assert neg_eidi.get_splicescore_flank(6) == (63170914-2-1, 63170914+20)
    assert neg_eidi.get_splicescore_flank(6, False) == (63170914-2, 63170914+20)

    assert pos_eidi.get_ss_flank(1,)  == (59704208-150-1, 59704208+150)
    assert pos_eidi.get_ss_flank(2, ups_w=100) == (59704624-100-1, 59704624+150)
    assert pos_eidi.get_ss_flank(3,dws_w=100,based0=False) == (59704750-150, 59704750+100)
    assert pos_eidi.get_ss_flank(4,ups_w=100, dws_w=100,based0=False) == (59705728-100, 59705728+100)
    assert pos_eidi.get_ss_flank(5,ups_w=120, dws_w=100,based0=False) == (59705756-120, 59705756+100)

    assert neg_eidi.get_ss_flank(1,)  == (63170914-150-1, 63170914+150)
    assert neg_eidi.get_ss_flank(2, ups_w=100) == (63172154-100-1, 63172154+150)
    assert neg_eidi.get_ss_flank(3,dws_w=100,based0=False) == (63172257-150, 63172257+100)
    assert neg_eidi.get_ss_flank(4,ups_w=100, dws_w=100,based0=False) == (63175637-100, 63175637+100)

    assert neg_eidi.get_ss_flank(1,strand_sp=True)  == (63176428-150-1, 63176428+150)
    assert neg_eidi.get_ss_flank(2, ups_w=100,strand_sp=True) == (63175688-150-1, 63175688+100)
    assert neg_eidi.get_ss_flank(3,dws_w=100,strand_sp=True,based0=False) == (63175637-100, 63175637+150)
    assert neg_eidi.get_ss_flank(4,ups_w=100, dws_w=110,strand_sp=True,based0=False) == (63172257-110, 63172257+100)


def test_RI_SuppaEventCoord():
    pos_eid = "ENSMUSG00000026127.13;RI:chr1:34443769:34443923-34443997:34444091:+"
    neg_eid = "ENSMUSG00000026121.13;RI:chr1:36551726:36551838-36551927:36552149:-"

    pos_eidi = SuppaEventCoord(pos_eid)
    neg_eidi = SuppaEventCoord(neg_eid)

    assert pos_eidi.get_alter(False)  == (34443923, 34443997)
    assert pos_eidi.get_alter()  == (34443923-1, 34443997)

    assert neg_eidi.get_alter(False) == (36551838, 36551927)
    assert neg_eidi.get_alter() == (36551838-1, 36551927)


    assert pos_eidi.get_ss_range(1, 4)  == (34443769-1, 34444091)
    assert pos_eidi.get_ss_range(2, 3, based0=False) == (34443923, 34443997)
    assert pos_eidi.get_ss_range(1, 3, True, False) == (34443769, 34443997)

    assert neg_eidi.get_ss_range(2, 4)  == (36551838-1, 36552149)
    assert neg_eidi.get_ss_range(2, 3, based0=False) == (36551838, 36551927)
    assert neg_eidi.get_ss_range(1, 3, True, False) == (36551838, 36552149)

    assert pos_eidi.get_splicescore_flank(1)  == (34443769-20-1, 34443769+2)
    assert pos_eidi.get_splicescore_flank(1, based0=False) == (34443769-20, 34443769+2)
    assert pos_eidi.get_splicescore_flank(2) == (34443923-2-1, 34443923+6)
    assert pos_eidi.get_splicescore_flank(2, False) == (34443923-2, 34443923+6)
    assert pos_eidi.get_splicescore_flank(3) == (34443997-20-1, 34443997+2)
    assert pos_eidi.get_splicescore_flank(3, False) == (34443997-20, 34443997+2)
    assert pos_eidi.get_splicescore_flank(4) == (34444091-2-1, 34444091+6)
    assert pos_eidi.get_splicescore_flank(4, False) == (34444091-2, 34444091+6) 

    assert neg_eidi.get_splicescore_flank(1)  == (36552149-2-1, 36552149+20)
    assert neg_eidi.get_splicescore_flank(1, based0=False) == (36552149-2, 36552149+20)
    assert neg_eidi.get_splicescore_flank(2) == (36551927-6-1, 36551927+2)
    assert neg_eidi.get_splicescore_flank(2, False) == (36551927-6, 36551927+2)
    assert neg_eidi.get_splicescore_flank(3) == (36551838-2-1, 36551838+20)
    assert neg_eidi.get_splicescore_flank(3, False) == (36551838-2, 36551838+20)
    assert neg_eidi.get_splicescore_flank(4) == (36551726-6-1, 36551726+2)
    assert neg_eidi.get_splicescore_flank(4, False) == (36551726-6, 36551726+2)


    assert pos_eidi.get_ss_flank(1,)  == (34443769-150-1, 34443769+150)
    assert pos_eidi.get_ss_flank(2, ups_w=100) == (34443923-100-1, 34443923+150)
    assert pos_eidi.get_ss_flank(3,dws_w=100,based0=False) == (34443997-150, 34443997+100)
    assert pos_eidi.get_ss_flank(4,ups_w=100, dws_w=100,based0=False) == (34444091-100, 34444091+100)

    assert neg_eidi.get_ss_flank(1,)  == (36551726-150-1, 36551726+150)
    assert neg_eidi.get_ss_flank(2, ups_w=100) == (36551838-100-1, 36551838+150)
    assert neg_eidi.get_ss_flank(3,dws_w=100,based0=False) == (36551927-150, 36551927+100)
    assert neg_eidi.get_ss_flank(4,ups_w=100, dws_w=100,based0=False) == (36552149-100, 36552149+100)

    assert neg_eidi.get_ss_flank(1,strand_sp=True)  == (36552149-150-1, 36552149+150)
    assert neg_eidi.get_ss_flank(2, ups_w=100,strand_sp=True) == (36551927-150-1, 36551927+100)
    assert neg_eidi.get_ss_flank(3,dws_w=100,strand_sp=True,based0=False) == (36551838-100, 36551838+150)
    assert neg_eidi.get_ss_flank(4,ups_w=100, dws_w=110,strand_sp=True,based0=False) == (36551726-110, 36551726+100)
