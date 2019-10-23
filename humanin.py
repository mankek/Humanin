import pandas as pd


def add_aa(dict_in, set_in):
    for place in dict_in.keys():
        for i in set_in:
            if i not in dict_in[place].keys():
                dict_in[place][i] = 0
    return dict_in


def count_positions(text_in, text_dict, set_in):
    for place, acid in enumerate(text_in):
        text_dict[place].setdefault(acid, 0)
        text_dict[place][acid] += 1
        set_in.add(acid)
    return text_dict, set_in


def read_proteins():
    aa_set = set()
    pos_dict = dict()
    for i in range(0, 24):  # 24: the number of amino acids in the peptide
        pos_dict[i] = dict()
    with open("humanin.txt", "r") as file_in:
        for i in file_in.readlines():
            pos_dict, aa_set = count_positions(i.rstrip(), pos_dict, aa_set)
    pos_dict = add_aa(pos_dict, aa_set)
    return pos_dict, aa_set


pos_dict_out, aa_set_out = read_proteins()
print(pos_dict_out)


def get_counts(pos_dict, aa_set):
    counts_dict = dict()
    for aa in aa_set:
        counts_dict.setdefault(aa, [])
    for place in pos_dict.keys():
        for aa in pos_dict[place].keys():
            counts_dict[aa].append(pos_dict[place][aa] + 1)
    counts_df = pd.DataFrame(counts_dict)
    probs_df = counts_df.applymap(lambda x: (x/32) if x != 0 else x)  # 32: twice the number of sequences considered to find a consensus
    consen_dict = {"Sequence": probs_df.idxmax(1), "Probabilities": probs_df.max(1)}
    consen_df = pd.DataFrame(consen_dict)
    consensus = "".join(counts_df.idxmax(1).values)
    return probs_df, consen_df, consensus


probability_df, consensus_df, consensus_text = get_counts(pos_dict_out, aa_set_out)


def find_prob(prob_df, text_in):
    aa_list = list(prob_df.columns)
    text_probs = 1
    for index, acid in enumerate(text_in):
        if acid in aa_list:
            prob = prob_df.at[index, acid]
            text_probs = text_probs * prob
        else:
            prob = 1/32
            text_probs = text_probs * prob
    return text_probs


def find_humanin(sequence_in):
    prob = 0
    likely_hum = ""
    sequence_in.replace("-", "*")
    for i in range(0, len(sequence_in) - 24):
        subseq = sequence_in[i:(i + 24)]
        sub_prob = find_prob(probability_df, subseq)
        if sub_prob > prob:
            prob = sub_prob
            likely_hum = subseq
    return prob, likely_hum


test_seq_1 = "TKNTPAQPLQPSQNTN-NITTH-VYEMELYLKRNSM-YLKGIMKESIYL-VYHS-TYTPYLLHHGLASMLQI-RPTPVTPKPDELPKSNRKDEPISVAMEWEDF-VEAKHQPSLVMAGYSLNEL-FNSSHTPNMPQFNGTA--YTMKVQLYCNWLQPKLVGHTFLNMTVGL-AATKNKCVQAHMTPMQ-FHHTPY-PNWVILLLNK-NPAKTSNKSPISSQALKSARTNHW-LTDYKKYYKQNQVYLPLTLLSPQH-HA--MMKSLK-NSANNTPTVYQKHSL-LQSSIKGPACPVKMLF-RPRYPNRAKVA-SIAP-MGAGMNGWM-VYLSPETN-WNWSFSTKAGIHP-DEQTLWSLKL-YYNQFTPTFPGNSLKLDSINTFSWGDCGM-KIFHVN-TMIHP-PTSQ-KTHLTQ-HWSTNQVTPGMTAPSPS-AHIAKGAYDLDVGSGQPSGAAATKGSFVQRLMSYVIWVQTGEIQVGFYLCLSFLSTKGPEKLTPLLTNTRLISLLNTNQNTKTILLTPLDKG--"
test_seq_2 = "PKMRLHNRSNHHKTQIKTLLHTKYM--NFILSAMAYSTS-ESWKNPSTFKYITAGLTPRTSCIMV-QAYY-FSALHLSPRNQTSYL-AIV-MNPSLWQ-SGKTS---RNTNRAWW-LVTH-MNFNSTLATPPMYHNLMAQLEDMQ--YSFIVIGYNQN-WDTLSSM-LWAFKPPPKTNAFKHM-HQYNNFTTPPMDQTELFYYSMKETLLKLVMSHPSPRKPLSQHGPTTDS-QTTKSIMNKIKYIYHLHSCHPNTGMHKG--KVSKGTRLMTPQLFTKNMAFSYNQVLKVPPAQWKYFFNGRGILTVQ--RNQLPLK-GLVWTAEWGCTCLL-PINETDLSVQKLEFTHKTN-PCGA-NSNTTTNLPLLSPETP-NLMVLMLLVGATAEYKKSST-TGP-SIQDQQVKEKHIWPSNTDQRTKLPQG-QRHPLQEPMSP-GLTTSMLDQDNQVAQPLLKVRLFND-CPTWSEF-PEKS-SVSIYV-VSSVRKDQ-N-RHY-PTRV-FHYWTQIKTPKQSF-HP-MKG-"
test_seq_3 = "QKYACTTAPTITKHKLKHYYTLSMWD-TLS-AQ-HMVPQGNHE-IHLPLSMSQQDLHPVPLASWFSKHTTDLAPYTCHPET-RAT-EQS-GWTHLCGN-VG-LLG-GETPTEPGDSWLLIKWTLIQL-PHPQYTTI-WHS-KMYNKGTALL-LATTKISGTHFPQYNCGPLSRHQKQMRSSTYNTNTMISPHPL-TKLSYFITQ-KKPC-N---VTHLLASP-VSTDQPLMVN-LQKVL-TKSSMFTTYTPVTPTQACMKDNKKSQKELG--HPNCLPKT-PLATIKY--SRLPSENTFLTAAVS-PCKGSVINCPLN-GWYERLNEGVPVSWDQLMKLIFQYKSWNSPM-RTDPVELKTLMLQPIYPYFPRKLLKTW-Y-YF-LGRLRNMKNLPRKQDHNPSKTNKSKKNTSDPVTLINEPSYP-DNSAIPFKSPYRQGGLRPRCWI-TTKWRSRY--FVCSTINVLRDLSSDR-NPGRFLSMFKFPQYE-T-ETNATINQHAFNFTTEHKSKHQNNPFNTP---VE"
test_seq_1_rev = "SLPFI-GC-KDCFGVLICVQ-WN-TRVG--WR-FLWSFRTEET-T-METDLDFSGLNSDHVGH-SLNKRTFSSGCATWLSWSNIEVVSPLGDMGSW-GWRCYPWGNLVRWSVLLGQMCFSLTCWSWMDYGPVYVEDFLYSAVAPTKSINTIKF-GVSGES-GKLVVVLEF-APQGLFVLWVNSSFCTE-SVSLIGL--QVHPHSAVHTSPYL-GNWLRYLCTV-MPRPLKKYFHWAGGTFNTWL-LKAMFLVNSWGVISRVPFETFYYPLCMPVLGWQECKW-MYLILFMMLFVVC-LSVVGPCWLKGLRGDGWLITSFS-VSFIE--NNSVWSMGGVVKLLYWCYMCLNAFVFGGGLKAHSYIEESVSH-FWL-PITMKLYLYCMSSSCAIKLWYIGGVA-VELKFI-WVTSYHQARLVFRLYLEVFPLYCH-DGFILTIAL--LVWFRGD-C-ALNL-YAC-TMMQEVRGVSPAVMYLKVDGFFHDSLEVLYAIALKMKFYLMYLVCSNVLICVLWWLERLC-RIFG"
test_seq_2_rev = "LYPLS-GVK-IVLVFWFVFSSEIKRVLVNSGVSFSGPFVL-KLKH--KPTWISPVWTQIT-DINRWTNEPLVAAAPLGCPDPTS-S-APLAMWALEGDGAVIPGVTWFVDQCYWV-CVFLWLVGLGWIMVLFTWKIFYIPQSPQLKVLMLSSFKEFPGKVGVNWL-Y-SFKLH-VCSSYGWIPAFVLKDQFH-LVSGD-YTLIQPFMPAPI-GAIDYATFARLGYRGR-KSIFTGQAGPLMLDCS--LCFW-TVGVLLAEFLL-LFIILYACLCWGD-SVSGKYTWFCL-YFL-SVNYQWLVRADL-ACEEMGDLLLVLAGFLLLSNKMTQFGL-GVWWNYCIGVMCAWTHLFLVAA--PTVML-KVCPTNFGCSQLQ-SCTFIVYLLAVPLNCGMLGVWLELN-SSFNE-PAIT-LGWCFAST-KSSHSIATEMGSSLRLLLGSSSGFGVTGVGR-ICSMLAKPWCK-YGV-VLLWYT---MDSFMIPL-YYMLLRL--SSISYT-CVVMF-FVFCDGWSGCAGVFL"
test_seq_3_rev = "STLYLGVLKGLFWCFDLCSVVKLNACWLMVALVSLVLSYWGNLNMD-NRPGFLRSEL-SR-TLIVEQTNL--RLRHLVVLIQHRGRKPPWRYGLLKGMALLSLG-LGSLISVTGSDVFFFDLLVLDGLWSCLRG-FFMFRSRPN-KY-YYQVL-SFRGK-G-IGCSI-VLSSTGSVRLMGEFQLLYWKISFINWSQETGTPSFSRSYQPLFKGQLITLPLHG-DTAAVKKVFSLG-RDL-YLIVAKGYVFGKQLGCY-PSSFWDFLLSFMHACVGVTGV-VVNMLDFVYNTFCSLLTISGWSVLT-GLA--WVTYY-F-QGFFYWVMK-LSLVY-GCGEIIVLVLYVLERICFWWRLKGPQLYWGKCVPLILVVANYNKAVPLLYIF-LCH-IVVYWGCG-SWIKVHLMSNQLSPGSVGVSPLP-SLPTLLPQ-WVHPYDCS-VARLVSGWQV-GAKSVVCLLNHDA-GTGCKSCCDMLKG-WILSWFPWGTMCYCA-DKVLSHMLSV--CFNLCFVMVGAVVQAYFW"

gracilis_1 = "A-SSPHPPPPPTTKTTKPTKQTIWQGQ-MRLNPTNDAMVKQYRKGKVKYPQL-SNMKQGQTPVPFASWSSKNPLAK-IYSPLPRNRMSYFWATN-PNPSL-QKSGTTPK-GQKT-RTRW-LVTW-KNISSTS-RTTH-LVLTWRL--YSMGVQLYWT-LQPAAENNQNHAISGPQSSHHQTTASQPHHKNTTLNTSPLQAHQAIL-PY--TYA-ISNKSPSLSAHACNSAMEHLLTIN-PKKDQCYTNKTKILPTVNPTQDCLM-KIKHQ--N-ANFQAPTVYQKHSL-QKVLKVSPAQWLTTYLNGRGILTVQ--RNHLSPK-GLVWMAKWGPSCLPWLISEIDLPVQKLVPCGA-NTYAKQHQPPNNAMTFSVGATTETNKTSKP-ISSSTTDKHVYITWP-TTLNNEPSYP-DNSAIFFKSPHRQEGLRPRCWI-TPQWCSRY--FVCSTIKVLRDLSSDR-NPGRFLSMLNLPQYE-T-KVGPILT-SP-STADMNLTCR-QKSKPKT-AI"

human_test_3 = "-T-PQTHSTLLPDNLSQTIYPNKV-AIEIETWRNRYSTARER-KIITKHNIARTNPYTFCIMN-LEITLQGEPKLRPPKPDELPKNS-KSTPVYVAK-WEDL-VEATNLPSLVIAGCPR-NLSSTLNLPTEPSKSPCKFNC-SKEEQLFGH-EKTL-RE-KI-HP--A-KQPPIKKAFKLNTHYLKNPKHITELLTPNWTNLSPYRRTNVSISNMKTFSSA-ACVRLKH-TDN-QPNIYNQPTSHYYPHCQPNTGMLIRKG-KK-KELGKSYPACLPKTSPLASPVLEAPPAQ-HMFNGRGTLTVQR-HNHLFLK-GPV-MAPRGFSCLLLLTSEIDLPVKRRA-HSKTRRPYGALIY-CKQYLTNPQVLNYQTCIKNFGWGDLGAEPNLRAVHAKTSPVKANYYTQLIQ-LDQRNKLP-G-QRNPILESISTIGFTTSMLDQDIPMVQPLLKVRLFND-SPT-SEFRPE-SRSVSIYXQIPPCTKGQEK-GLLHKAPSPVNDIIST-YYTHTHPRTGF"

print(find_prob(probability_df, consensus_text))
print(consensus_text)
print("Gecko frame 1: ", find_humanin(test_seq_1))
print("Gecko frame 2: ", find_humanin(test_seq_2))
print("Gecko frame 3: ", find_humanin(test_seq_3))
print("Gecko reverse frame 1: ", find_humanin(test_seq_1_rev))
print("Gecko reverse frame 2: ", find_humanin(test_seq_2_rev))
print("Gecko reverse frame 3: ", find_humanin(test_seq_3_rev))
print("Glass lizard frame 1: ", find_humanin(gracilis_1))
print("Human frame 3: ", find_humanin(human_test_3))


