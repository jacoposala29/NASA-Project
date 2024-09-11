# jupyter notebooks 
notebooks for creating composite plots for MHW events: number of events, duration, intensity, heat budget terms composites, first/second/third contributor, onset/decline contributors. Line plots for horizontal average in regions


these are all the run needed for composite plots (made by ECCOv4r4_load_mat_files_and_create_maps_DAILY.ipynb):

########### OISST MONTHLY - DONE
1993-2016 OISST monthly SST - DONE
case_sel = get_dict_4_case_of_interest(tag_case = 'OISST MONTHLY 1993-2016 plot years 1993-2016')

2004-2016 OISST monthly SST - DONE
case_sel = get_dict_4_case_of_interest(tag_case = 'OISST MONTHLY 2004-2016 plot years 2004-2016')

1992-2017 OISST monthly SST - DONE
case_sel = get_dict_4_case_of_interest(tag_case = 'OISST MONTHLY 1992-2017 plot years 1992-2017')


########### ECCO MONTHLY - DONE
1993-2016 ECCO monthly k0-k5 - DONE
case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO MONTHLY OHC k0-k5 1993-2016 plot years 1993-2016')

2004-2016 ECCO monthly k0-k5 - DONE
case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO MONTHLY OHC k0-k5 2004-2016 plot years 2004-2016')

1992-2017 ECCO monthly k0-k5 - DONE
case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO MONTHLY OHC k0-k5 1992-2017 plot years 1992-2017')

1993-2016 ECCO monthly zlev01 - DONE
case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO MONTHLY zlev01 1993-2016 plot years 1993-2016')

2004-2016 ECCO monthly zlev01 - DONE
case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO MONTHLY zlev01 2004-2016 plot years 2004-2016')

1992-2017 ECCO monthly zlev01 - DONE
case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO MONTHLY zlev01 1992-2017 plot years 1992-2017')


########### ARGO MONTHLY - DONE
2004-2016 Argo monthly k0-k5 - DONE
case_sel = get_dict_4_case_of_interest(tag_case = 'ARGO MONTHLY OHC 2004-2016 plot years 2004-2016')

2004-2017 Argo monthly k0-k5 - DONE
case_sel = get_dict_4_case_of_interest(tag_case = 'ARGO MONTHLY OHC 2004-2017 plot years 2004-2017')


########### OISST - DAILY gap 0 and gap 2 - DONE
 1992-2017 OISST daily k0 (5+) gap 0 - DONE 
 case_sel = get_dict_4_case_of_interest(tag_case = 'OISST DAILY k0 1992-2017 5+days plot years 1992-2017 gap 0')

 1992-2017 OISST daily k0 (5-30) gap 0 - DONE 
 case_sel = get_dict_4_case_of_interest(tag_case = 'OISST DAILY k0 1992-2017 5-30days plot years 1992-2017 gap 0')

 1992-2017 OISST daily k0 (30+) gap 0 - DONE 
 case_sel = get_dict_4_case_of_interest(tag_case = 'OISST DAILY k0 1992-2017 30+days plot years 1992-2017 gap 0')

 1992-2017 OISST daily k0 (5+) gap 2 -  DONE
 case_sel = get_dict_4_case_of_interest(tag_case = 'OISST DAILY k0 1992-2017 5+days plot years 1992-2017 gap 2')

 1992-2017 OISST daily k0 (5-30) gap 2 -  DONE
 case_sel = get_dict_4_case_of_interest(tag_case = 'OISST DAILY k0 1992-2017 5-30days plot years 1992-2017 gap 2')

 1992-2017 OISST daily k0 (30+) gap 2 -  DONE
 case_sel = get_dict_4_case_of_interest(tag_case = 'OISST DAILY k0 1992-2017 30+days plot years 1992-2017 gap 2')


########### ECCO DAILY k0 - gap 0 and gap 2 - DONE
 1992-2017 ECCO daily k0 (5+) gap 0 -  DONE
 case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO DAILY zlev00 1992-2017 5+days plot years 1992-2017 gap 0')

 1992-2017 ECCO daily k0 (5-30) gap 0 -  DONE
 case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO DAILY zlev00 1992-2017 5-30days plot years 1992-2017 gap 0')

 1992-2017 ECCO daily k0 (30+) gap 0 -  DONE
 case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO DAILY zlev00 1992-2017 30+days plot years 1992-2017 gap 0')

 1992-2017 ECCO daily k0 (5+) gap 2 -  DONE
 case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO DAILY zlev00 1992-2017 5+days plot years 1992-2017 gap 2')

 1992-2017 ECCO daily k0 (5-30) gap 2 -  DONE
 case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO DAILY zlev00 1992-2017 5-30days plot years 1992-2017 gap 2')

 1992-2017 ECCO daily k0 (30+) gap 2 -  DONE
 case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO DAILY zlev00 1992-2017 30+days plot years 1992-2017 gap 2')


########### ECCO DAILY OHC k0-k5 - gap 0 and gap 2 - 
 1992-2017 ECCO daily OHC k0-k5 (5+) gap 0 -  DONE
 case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO DAILY OHC k0-k5 1992-2017 5+days plot years 1992-2017 gap 0')

 1992-2017 ECCO daily OHC k0-k5 (5-30) gap 0 -  DONE
 case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO DAILY OHC k0-k5 1992-2017 5-30days plot years 1992-2017 gap 0')

 1992-2017 ECCO daily OHC k0-k5 (30+) gap 0 -  DONE
 case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO DAILY OHC k0-k5 1992-2017 30+days plot years 1992-2017 gap 0')

 1992-2017 ECCO daily OHC k0-k5 (5+) gap 2 -  DONE
 case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO DAILY OHC k0-k5 1992-2017 5+days plot years 1992-2017 gap 2')

 1992-2017 ECCO daily OHC k0-k5 (5-30) gap 2 - DONE 
 case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO DAILY OHC k0-k5 1992-2017 5-30days plot years 1992-2017 gap 2')
 
 1992-2017 ECCO daily OHC k0-k5 (30+) gap 2 -  DONE
 case_sel = get_dict_4_case_of_interest(tag_case = 'ECCO DAILY OHC k0-k5 1992-2017 30+days plot years 1992-2017 gap 2')

