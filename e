[33mcommit f1a7a86cf86daf3ec92d04630ff4550e3a2fec65[m[33m ([m[1;36mHEAD -> [m[1;32mFRED-covid19-usa[m[33m, [m[1;31morigin/HEAD[m[33m, [m[1;31morigin/FRED-covid19-usa[m[33m)[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Sat May 22 16:24:05 2021 -0400

    Fixing issue with priority phases. Do not shuffle after adding to overall queue when priorities are enabled

[33mcommit 34c834b7cee1fc3efc0df029bb7122ca220daa53[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed May 19 22:56:17 2021 -0400

    Bypassing efficacy against infection when checking for vaccines availabel for a person in the priority queue

[33mcommit 4cee31338f0b92b00adc96f84949be739d9e8d20[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Mon May 17 13:24:03 2021 -0400

    Increasing speed of vaccine manager removal of dead people

[33mcommit 506b6221b1e10992445f12583202d44ab10fb045[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Sat May 8 00:33:02 2021 -0400

    Improving vaccination performance

[33mcommit 5161c1e39e99ad3a80ac8897b871d0f32f8294f8[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed May 5 21:23:59 2021 -0400

    Added a lap for vaccine update

[33mcommit 272483a3c8e0539c01f4a5297830ab45c6f19d9f[m
Author: Guido <confunguido@gmail.com>
Date:   Tue May 4 21:54:21 2021 -0400

    Remove unnecessary prints and fix max parameters issue

[33mcommit 3cbdd864fba11403a59d09d264225a560a247bd9[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Tue May 4 19:45:33 2021 -0400

    Allow cross protection for vaccines

[33mcommit 211bbe7ae0ddac4760becb9531717aef726f4674[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Sun May 2 21:40:52 2021 -0400

    Nursing home added to place when imported

[33mcommit d9d3348ce7233489142658ac75f8be4378582d23[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Sun Apr 25 16:07:03 2021 -0400

    Hospitalizations can be multiplied by a scaling factor that varies in time (if enabled). Keeping track of total number of infections an agent has

[33mcommit 0a1c3bb274d2bfb207a28a6d038e2f69ba47641e[m
Author: Guido <confunguido@gmail.com>
Date:   Tue Apr 20 21:17:31 2021 -0400

    Fixing issue with immunity waning in the presence of various diseases with cross-protection

[33mcommit 40d6dc69a32c77640c0852c0732a4044ca26b721[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Tue Apr 13 21:02:47 2021 -0400

    Add immunity loss to only a proportion of the population

[33mcommit c2d933f50625a8797c3a5551df688700096a60ac[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Fri Apr 9 10:18:26 2021 -0400

    fix immunity duration

[33mcommit 434b5f5dcf6e63b0feacc1cbec54a000b4f73a16[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Tue Apr 6 12:57:42 2021 -0400

    Minor adjustments for school reopening

[33mcommit cea93730765107448c65a3336a2b40b6425e7f64[m
Author: Guido <confunguido@gmail.com>
Date:   Tue Apr 6 10:51:29 2021 -0400

    Trying to fix School bug

[33mcommit 813626882fd1ef021f3c7fcb198bfe455fa2ab8d[m
Author: Guido <confunguido@gmail.com>
Date:   Mon Apr 5 09:58:32 2021 -0400

    School reopening priority to grade reopening

[33mcommit cd5f51b9265802ff1d4ca83fd9d93433acb4ddfd[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed Mar 31 15:46:30 2021 -0400

    Stop test aborting accidentally on

[33mcommit 597ced7c8a9efaa300a4097b70d501dc22ef4f39[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Fri Mar 19 19:57:50 2021 -0400

    Fixing hospitalization periods

[33mcommit 303f117a087d785e0173e264bce53c30224da9ab[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed Mar 17 12:39:17 2021 -0400

    Incorporating vaccination rollout by phases with different prioritization

[33mcommit 278b5b5b867ad60d5c5af15ec4efdba5b7b0a0ff[m
Merge: bd75c7b d180be6
Author: Guido Espana <confunguido@gmail.com>
Date:   Sun Mar 14 15:59:55 2021 -0400

    Merge branch 'FRED-covid19-usa' of github.com:confunguido/FRED into FRED-covid19-usa

[33mcommit bd75c7bd59524eb79c42bfec8b61d14fef8dfd53[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Sun Mar 14 11:48:27 2021 -0400

    Enable timeseries of community contact scaling factor

[33mcommit bf211a718c932f0a2e4f39dc7aa998374cc83562[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Sat Mar 13 16:36:02 2021 -0500

    Implementing vaccine capacity variable by day. Hospitalizations. Also included VE for symptoms and severe disease

[33mcommit d180be6401cd11048931105c66eb4a0766296a97[m
Author: Guido <confunguido@gmail.com>
Date:   Tue Feb 16 18:05:12 2021 -0500

    Adding school labels to infection output

[33mcommit 94b4a362ad221c02717cdf81e2277234b08e548b[m
Author: Guido <confunguido@gmail.com>
Date:   Thu Jan 21 14:46:48 2021 -0500

    Add classroom in infection.cc output fro income and classroom label

[33mcommit df97511e5af2a5de83d99d19ff08c63d6e575a7a[m
Merge: 6de4a4d 70f2779
Author: Guido <confunguido@gmail.com>
Date:   Wed Jan 20 11:10:39 2021 -0500

    Fixing merging issues

[33mcommit 6de4a4df5cadb2c902a16d02644f90a33cccc21e[m
Author: Guido <confunguido@gmail.com>
Date:   Wed Jan 20 11:06:13 2021 -0500

    School census tract in output of infections

[33mcommit 70f27798983d4e9d39d25a0045330fe3f825cff1[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Tue Jan 19 20:10:10 2021 -0500

    Add classroom label to infection output

[33mcommit 343157b45016e4159938ee31b926dacffb35e6c3[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Mon Jan 4 21:00:08 2021 -0500

    Adding holiday (diff. than seasonal) increase in contacts within a neighborhood

[33mcommit ff0c693f6f5c1290ea921fb2ee8266dc48afe434[m
Author: Guido <confunguido@gmail.com>
Date:   Mon Jan 4 13:22:22 2021 -0500

    Mixing group label in output

[33mcommit 566aa21ed0121d806d76e583e8ed9e62ed909384[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Mon Dec 21 15:50:40 2020 -0500

    Adding school income to infections output

[33mcommit 522d4adcb8e03527e66d2a2fd829b5380ad72e21[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Tue Nov 24 15:32:00 2020 -0500

    Fixed census tract id of schools

[33mcommit 4c9f1b6058d9f922be342fef67f794f9169e09e6[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Tue Nov 17 20:27:23 2020 -0500

    adding the possibility to close schools by census tract

[33mcommit 3f52d5536fdc2ac6f6608adc570abcba45aa2829[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Mon Oct 19 13:50:15 2020 -0400

    Adding classroom size by grade

[33mcommit 51db70f5f1521d47415c0db04e68585282d8f980[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed Oct 14 17:54:04 2020 -0400

    Nursing home importations have an infector now

[33mcommit b5ebbd8041250c188ca7a2b72737da2f7411245e[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed Oct 14 15:59:39 2020 -0400

    Adding mixing group info for nursing home importations

[33mcommit 5d00863d14685db5f14c4c7107524a86426fc851[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Thu Oct 1 19:17:04 2020 -0400

    Minor changes. Allowing classrooms for grades of people older than GRADES

[33mcommit f141843ed55cab1fd684816dc2170a31cc0ce886[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed Sep 30 23:41:42 2020 -0400

    Fixing issue with schools not being read

[33mcommit 1d89feb7b61a302e4128707748eb192ab9d2db0d[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed Sep 30 21:11:30 2020 -0400

    Adding the possibility for schools to have an income in the input

[33mcommit 4648c04a8a7c9b25d45089149d09fcfa7bc59012[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed Sep 30 17:41:48 2020 -0400

    Fixing an issue when reading capacity for schools

[33mcommit bee52669aa4b5f64337f7563a90707aec8f7fb45[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed Sep 30 14:55:09 2020 -0400

    Schools can be closed with a new global_schedule option, specifying start end grade_min grade_max capacity

[33mcommit 62cddf7a92cc88199054fa04eb266129d898927c[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Sun Sep 27 16:15:12 2020 -0400

    Fixed student/teacher ratio when gq are disabled. Added minimum age to wear a face mask

[33mcommit a6075015f71cc36aa276da0facc2d9137c263d07[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Tue Sep 22 22:39:54 2020 -0400

    Fixed bug on sheltering in place... full city was overriding census tract compliance

[33mcommit 19aaa5d8f144740804ab44ecf232271bd8d5bdd6[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Tue Sep 22 14:47:21 2020 -0400

    Default shelter in place file

[33mcommit 9feddc59c099254db16302266517747bb94aed12[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Mon Sep 21 23:25:32 2020 -0400

    Re-structure of shelter in place when using a file. Now houses can be sheltered for each census tract. If census tract is ingored, then -1 is used

[33mcommit a7f7d7bf4c2a3151b01fcabf6c14e6ca0f9928ba[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed Sep 9 10:31:06 2020 -0400

    Adding differential susceptibility threshold in addition to logistic susceptibility

[33mcommit 196a16c820bfd540b2c4eab9aa5c503e290106b0[m
Merge: af5c20e 2ef1a32
Author: Guido Espana <confunguido@gmail.com>
Date:   Mon Sep 7 15:37:00 2020 -0400

    Merge branch 'FRED-covid19-usa' of github.com:confunguido/FRED into FRED-covid19-usa

[33mcommit af5c20e220cbed0d61d3fd918febff847b5f6799[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Mon Sep 7 15:35:31 2020 -0400

    Fixed an issue preventing gnuplot to create the image

[33mcommit 2ef1a3288009b8af5d3831f800bfb3bc9861b788[m
Author: Guido <confunguido@gmail.com>
Date:   Thu Aug 13 11:35:35 2020 -0400

    Export id of house in infection list

[33mcommit 889dfbfc627559f6b5c5f89ba023cccf0418c177[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Mon Aug 10 23:03:04 2020 -0400

    Adding a student flag to the infections output

[33mcommit e925ae7db00f98884e9af895ba932064ca136323[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Mon Aug 10 17:14:57 2020 -0400

    Adding a high value for susceptibility by age

[33mcommit 52b655c30b6a28c8e082f7f1442250b028026f8c[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Fri Aug 7 15:26:11 2020 -0400

    Fixing teachers and students assignment

[33mcommit 5b3450909e666930682a375e1d9aa3a94f460b28[m
Author: Guido <confunguido@gmail.com>
Date:   Thu Jul 30 17:12:26 2020 -0400

    Fixing reduced school capacity

[33mcommit 6dca9a6f3854b8beba476edc90db3fdbac8a6d5c[m
Merge: 619187c d8d78db
Author: Guido <confunguido@gmail.com>
Date:   Thu Jul 23 15:30:01 2020 -0400

    Merging facemask branch

[33mcommit d8d78dbe82da8fbb0541f8d4faa71a52d7d7e806[m[33m ([m[1;31morigin/FRED-covid19-usa-face-masks[m[33m)[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Thu Jul 23 15:19:30 2020 -0400

    Adding prevalence of infected but not symptomatic

[33mcommit 99641d81f84bf1a8f2482dacad5165dd4c7a2236[m
Merge: d3479d1 8f71bf2
Author: Sean Cavany <scavany@nd.edu>
Date:   Thu Jul 23 12:04:09 2020 -0400

    Merge branch 'FRED-covid19-usa-face-masks' of https://github.com/confunguido/FRED into FRED-covid19-usa-face-masks

[33mcommit d3479d1bc79a515025a8c2d389391007efeb5bed[m
Author: Sean Cavany <scavany@nd.edu>
Date:   Thu Jul 23 12:03:07 2020 -0400

    fixed tiny but important bug with log-normal implementation.

[33mcommit 8f71bf2bf4e185d3bfeab301fb202b1fcc643191[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed Jul 22 18:47:17 2020 -0400

    Adding school reduced capacity

[33mcommit fb7d0821b66ed289460d5bb9e0d33f9f860639fb[m
Author: Sean Cavany <scavany@nd.edu>
Date:   Wed Jul 22 10:27:08 2020 -0400

    Stopped isolating if no longer symptomatic

[33mcommit 3abd597d9411ffdb67b29d8cce3610edda810eeb[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Tue Jul 21 14:33:24 2020 -0400

    adding facemask compliance with a timeseries instead of fixed compliance

[33mcommit 48e2f9681f8cc9ccf660c7e6355c2266c82f1f0b[m
Author: Sean Cavany <scavany@nd.edu>
Date:   Tue Jul 14 21:05:37 2020 -0400

    Fixed bugs in face mask and hand washing

[33mcommit cad6e9e15f59d74262db7b585d335c074981d3da[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Mon Jul 13 16:59:19 2020 -0400

    Fixing some issues reading parameters for face masks

[33mcommit 3e394cdbe6973797bcab94a77e210ef7560b9476[m
Author: Sean Cavany <scavany@nd.edu>
Date:   Mon Jul 13 14:33:45 2020 -0400

    updated face mask functionality to allow location-specific compliance

[33mcommit 619187c3a70fb971ad7bd83045b870a05f3aabda[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Thu Jul 9 13:20:09 2020 -0400

    fixing some bugs with moving average

[33mcommit 1f3503f39579e2a213de7cdd4db7385d2b21e749[m
Author: Sean Cavany <scavany@nd.edu>
Date:   Thu Jul 9 10:37:38 2020 -0400

    Updated peak calculations for MMODS

[33mcommit 702983e71f9b5a4408ae8c75d7144ef8093d34ff[m
Author: Sean Cavany <scavany@nd.edu>
Date:   Thu Jul 2 16:46:59 2020 -0400

    Corrected double initialization of parameter

[33mcommit 02bdd3931c4324133d8464bfeeb40518b351c49e[m
Merge: 5240b6c 10fb06d
Author: Sean Cavany <scavany@nd.edu>
Date:   Thu Jul 2 16:38:53 2020 -0400

    Merge branch 'FRED-covid19-usa-adaptive-shelter' into FRED-covid19-usa-face-masks

[33mcommit 10fb06d6f3bd21af47f2df8c9357ee4fe18c1a36[m
Author: Sean Cavany <scavany@nd.edu>
Date:   Thu Jul 2 16:37:11 2020 -0400

    Updates to face mask functionality

[33mcommit 5240b6cf327f743a8075042cafb4fe1c6233fa21[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Mon Jun 22 22:40:27 2020 -0400

    Adding importations factor for nursing homes in defaults

[33mcommit 2f72bc5ecc286caa925469327b1d02d551a3bece[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Mon Jun 22 22:34:05 2020 -0400

    Force importations on nursing homes driven by the prevalence of infections

[33mcommit ddfdf87c3668af11273371c19403afabcfa92d7c[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed Jun 17 16:51:54 2020 -0400

    Change age susceptibility to follow a logistic function

[33mcommit 8cd51f7e052edc8fe0d67e9e53e7d363754374a1[m
Merge: 54995e3 b617ec1
Author: Guido Espana <confunguido@gmail.com>
Date:   Fri Jun 12 14:43:05 2020 -0400

    Merge remote-tracking branch 'origin/FRED-covid19-usa-adaptive-shelter' into FRED-covid19-usa

[33mcommit b617ec1ace12d21031cb506654c8e0e5c38ef901[m
Author: Sean Cavany <scavany@nd.edu>
Date:   Fri Jun 12 13:44:53 2020 -0400

    corrected small error in shelter threshold

[33mcommit a291b46ac94f619f6566a2a6c8efe4643bc20293[m
Author: Sean Cavany <scavany@nd.edu>
Date:   Fri Jun 12 13:07:31 2020 -0400

    corrected small error in shelter threshold

[33mcommit 54995e309e410c382a4d790ccc110afb6e16eb48[m
Merge: 5bf7099 ddd6e73
Author: Guido Espana <confunguido@gmail.com>
Date:   Fri Jun 12 12:02:07 2020 -0400

    Merging addaptive sheltering interventions with main covid branch

[33mcommit ddd6e73f94ef961c188eddb2c3510c660ff05b9e[m
Author: Sean Cavany <scavany@nd.edu>
Date:   Fri Jun 12 09:44:49 2020 -0400

    error in Global.cc

[33mcommit 5bf7099cb58c0c7c17ceeaeb4e2f663f9af5566c[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Thu Jun 4 12:21:55 2020 -0400

    Adding susceptibility by age with an exponential function (monomolecular)

[33mcommit c30bb84a71861bcd3a380c41699f667fb8560d5d[m
Author: Sean Cavany <scavany@nd.edu>
Date:   Wed Jun 3 17:38:18 2020 -0400

    enabled adaptive end to shelter options

[33mcommit bdd4e534d5e7b36ea0f27634d355c5037ad99d94[m
Author: Guido <confunguido@gmail.com>
Date:   Tue Jun 2 17:02:48 2020 -0400

    Adding population by age to the age output level 4

[33mcommit 244b1eade64aa3517f0ceb63e3deb9bef91ffd10[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed May 27 15:51:36 2020 -0400

    Fixing bugs in Epidemic.cc: Symptomatic people were not being added to a list required for some age-specific symptomatic reports

[33mcommit 3b91530e926a51a4bbe0f65baca634faa6b2ac3c[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed May 27 12:54:04 2020 -0400

    Adding symptomatic report by age

[33mcommit ac5275a987d91f1caec7ba344b10fa87de3e157f[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Thu Apr 30 23:19:25 2020 -0400

    Fixing shelter-in-place timeseries so that houses don't have to re-enroll every day

[33mcommit 2e5791f3a9bea7a8ec642a180705ea6f6259db28[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Thu Apr 30 20:10:34 2020 -0400

    Adding the possibility of specifying shelter-in-place orders as a time series

[33mcommit 4e4a9fc2efbb20bb3b25c51d80ecca919c829fd7[m
Merge: 0cc7fa3 a1f1be3
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed Apr 29 12:47:48 2020 -0400

    Merge branch 'FRED-covid19-usa' of github.com:confunguido/FRED into FRED-covid19-usa

[33mcommit 0cc7fa37ec18d97f4e6e3e868d2adc0320cc2fdc[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Wed Apr 29 12:47:39 2020 -0400

    Forcing fred to use an ID

[33mcommit a1f1be36bd444bf3a67d756b761f76c69fa966ae[m
Author: confunguido <confunguido@gmail.com>
Date:   Sun Apr 19 11:10:42 2020 -0400

    Commiting before pulling changes from another branch

[33mcommit 2b9bb54997e0df5cdf2421b2f323e756ebd7d853[m[33m ([m[1;31morigin/FRED-v2.12.0[m[33m)[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Fri Apr 10 20:53:00 2020 -0400

    Add the possibility of doing linear compliance for shelter in place

[33mcommit 0d6d3480a2f169f950722d54235471e37244b4e9[m
Author: Guido <confunguido@gmail.com>
Date:   Mon Apr 6 19:15:39 2020 -0400

    Report invidiual case fatality, if enabled

[33mcommit 0d62dfd7cfe572c6807174359401c7108ca5b38c[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Fri Mar 27 21:05:13 2020 -0400

    Adding stepwise shelter-in-place strategy

[33mcommit a44b1a755b38ed0dceca745598602fbfbd1fc9c1[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Fri Mar 27 11:28:02 2020 -0400

    Adding a floor for compliance different than 0 after period is over

[33mcommit 4d8caea53a5a7f387d8a95070e5366bf88fe46b9[m
Author: Guido Espana <confunguido@gmail.com>
Date:   Thu Mar 26 20:05:04 2020 -0400

    Adding shelter by age and with or without students

[33mcommit 2815c34e99511fac6c834a56efc2a03d3f35f060[m
Author: Guido <confunguido@gmail.com>
Date:   Wed Mar 25 16:48:54 2020 -0400

    Add the ability to set the duration of infectiousness as a CDF, even when infectious period set to offset symptoms

[33mcommit 8fda4a334f394abfc337fd4506a888e0d4a1bbaa[m
Author: Guido <confunguido@gmail.com>
Date:   Sun Mar 22 17:29:31 2020 -0400

    Allowing the epidemic to use epidemic_offset so that infections won't spread before that value

[33mcommit b8c85b078f1ecfdac0930b69914a14d7787a86ce[m
Author: Guido <confunguido@gmail.com>
Date:   Wed Mar 11 12:43:14 2020 -0400

    Fixing fred string size so msa with many counties can fit

[33mcommit 5f54d5213a02d01bbf2036dcc4af25a01b6ffcfa[m
Author: Guido <confunguido@gmail.com>
Date:   Sat Mar 7 15:53:59 2020 -0500

    Adding reporting of deaths by nusing homes

[33mcommit b04130d3c459e830a89846dfc4602fd25be0dafc[m
Author: Guido <confunguido@gmail.com>
Date:   Tue Mar 3 16:16:27 2020 -0500

    Allowing asymptomatics to recover

[33mcommit 75e40daf0a81ae66c7410d48f60e5170a31f788e[m
Author: Guido <confunguido@gmail.com>
Date:   Mon Mar 2 19:32:56 2020 -0500

    Fixing the issue with dead people becoming non_infectious. Now people are terminated before their disease gets cleaned

[33mcommit 365907b70fed956e40fcabedff707a9a716d0605[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Wed Nov 16 10:15:31 2016 -0500

    Fixed a bug for Sick Leave
    
    Sick Leave quick fix

[33mcommit 6a2d60cf837e7a46c69e2bd2816e95240ddf2166[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Jan 16 10:40:15 2016 -0500

    Memory reduction in Markov_Epidemic

[33mcommit 6d54027bb02f52e0c4dc5c0f81d9827f8b643698[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jan 15 17:41:31 2016 -0500

    Markov transitions updated when changing age group.

[33mcommit 7323733403937a89b3a4334c2a842bb8219badd1[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jan 15 17:20:45 2016 -0500

    Added interaction between markov models and population dynamics

[33mcommit 2695f28381c976a836df8aafdc0c5d261ff66c44[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jan 15 11:13:28 2016 -0500

    Reformatted markov model parameters

[33mcommit 6637756e8328691051e7c002a30555f197a4fbcf[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jan 15 11:02:35 2016 -0500

    Age-specific Markov Models
    
    Markov models now have age-specific initial state distributions and
    transition probabilities.

[33mcommit 60a752420afcb43435f815f1bbec0212653035e9[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jan 15 09:36:31 2016 -0500

    Visualization code tweaks

[33mcommit d1294a6b32ffe801f8b7ee8e6e3b967416a04699[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Jan 14 16:28:56 2016 -0500

    visualization code cleanup

[33mcommit a92a5c3ae73101b4032bf4636a016cdcbf93efe0[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Jan 14 11:15:12 2016 -0500

    Added bin/README.txt

[33mcommit 8456ff81d64f226a5636d17a303e6c7741c79573[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jan 13 17:28:16 2016 -0500

    Code cleanup

[33mcommit b58021288ee142f1442e9c568fa97eee642b3f81[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jan 13 17:10:05 2016 -0500

    Updated fred_get_googlemap
    
    Eliminates FRED dependency on R,

[33mcommit 37750d12770c79720586678d8a7b6bdcbd5264df[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jan 12 17:38:28 2016 -0500

    Visualization upgrades

[33mcommit 968b7adf6b98e86d21e0bffb5776e9e8e6cb355a[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jan 12 16:48:13 2016 -0500

    Visualization support for case_fatalities

[33mcommit d072a54d87aa17f2008f53895a59bee24b9fd178[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jan 12 12:45:45 2016 -0500

    Add total case fatalities TD to output.
    
    No change to rt.

[33mcommit ceba2d78e50ada0103ade0f70b5b5d019b4afc4d[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jan 12 11:27:31 2016 -0500

    BUG FIX: Infections terminated properly when agent dies.
    
    No change to rt.

[33mcommit 787b7794684a71f928a42e70667f91267bc77c7e[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Mon Jan 11 13:48:25 2016 -0500

    Incremental Ages & Stages
    
    Check in of some fixes / coding standardization for Age & Stages. The
    rt works, but more work is needed on long-term Dynamic Population.

[33mcommit f378dd75ec66ddbc01cb9b71b70716f3dafea400[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sun Jan 10 20:40:01 2016 -0500

    markov_epidemic update

[33mcommit bfd47fe5971cf0c6e88fcc9c34cddc78bbfd717b[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sun Jan 10 14:46:20 2016 -0500

    Updated regression tests after last bug fix.

[33mcommit 85171550476e98fb6672af6e78a559ec96b61b42[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jan 8 15:36:56 2016 -0500

    minor code cleanup

[33mcommit 27ead631bc8c6d668e73da70f9b734db99238d14[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jan 8 14:57:32 2016 -0500

    fixed recovery bug
    
    A person only recovers when both the symptoms period and the infectious
    period have expired.

[33mcommit f52f85e567d59fb680c24fd848beed36de9c506e[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jan 8 14:05:20 2016 -0500

    updates to markov models

[33mcommit dcc0872cb1f355bdd3d36c5e52e7f3a5c7f508f7[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Jan 7 16:34:21 2016 -0500

    reduced required parameters for markov natural history

[33mcommit a45e04ad99c865724a3c2b1a2d3fd979b3c6be88[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Jan 7 11:32:02 2016 -0500

    working markov epidemic model

[33mcommit 9bc2dce0ce5d5d3ee046f39208ceb614172a86ac[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jan 6 16:16:56 2016 -0500

    bug fixes to fred_job

[33mcommit 99f1e301945c37c475e6cfc5281d76d08302c8f4[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jan 6 12:36:16 2016 -0500

    fred_job supports recursive includes

[33mcommit 1f840a8031ef6ecc2a025587d85f4a2747df2001[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jan 6 11:37:08 2016 -0500

    fred_job now handles all file parameters correctly

[33mcommit 6a125d7fe40143e3407b6e806f8971444fb5376b[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Wed Jan 6 09:18:37 2016 -0500

    County Demographics Reporting
    
    Initial Check-in of County Demographics Reporting

[33mcommit 1ab8a8f20322a6cc2c954cdf4d7e0a9eb716f9f8[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Jan 5 16:37:35 2016 -0500

    support for markov disease process

[33mcommit 820574b7f37d31c410a7b3481ff10a705cb59510[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Jan 5 15:17:02 2016 -0500

    Support for Markov models

[33mcommit 64875d71c717cfd3220f8cac068218e5736de3cf[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Jan 5 14:54:28 2016 -0500

    Initial support for Markov health processes

[33mcommit 92a31fb441b6d1a56dad1f1a2b0e21c2c0301a8e[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Mon Dec 21 21:26:50 2015 -0500

    Updated support for Markov Epidemics

[33mcommit aee8009d752a44ce5547bc2409c6551e54a15813[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Dec 18 09:33:41 2015 -0500

    Support for Markov Infections

[33mcommit 57a18ccc7ca1181ecb1e62fc3bd6e158821438cd[m
Merge: c32dcb8 c6b9514
Author: David Galloway <gallowdd@hotmail.com>
Date:   Thu Dec 17 14:26:51 2015 -0500

    Merge pull request #40 from FRED-PITT/synch_sl_county
    
    Synch up Sick Leave and County demo

[33mcommit c6b95142b005b317713337767a1d1b781b3a358f[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Thu Dec 17 14:24:56 2015 -0500

    Synch up Sick Leave and County demo
    
    Just want to synch up my branch so I am not too far from main branch …
    Changed the way sick leave is distributed to workplaces
    Added some code for counties to track demographics
    Removed unnecessary Self pointers in Activity (it already has a pointer
    to the agent)

[33mcommit c32dcb8c339326b9e802b98268aeb2e9dee1d92d[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Dec 16 15:12:03 2015 -0500

    fixed typo

[33mcommit 2adfa29a9bba0eaf7c9d8f73a02519499cc7aa0a[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Dec 16 15:11:07 2015 -0500

    support for overdose deaths

[33mcommit 4556ce88329f976ce0149dc7cd6db1bca1f4c260[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Dec 16 14:57:56 2015 -0500

    Fixed initialization.

[33mcommit 5936f54f4f52b93f86b00a3bde7925199b5ff346[m
Merge: b213cd5 816a888
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Dec 15 17:24:06 2015 -0500

    Merge pull request #39 from FRED-PITT/visualization_update
    
    updated visualization counters

[33mcommit 816a888b07bb8a645a0413a2f813f1aef39d3f34[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Dec 15 17:23:17 2015 -0500

    updated visualization counters

[33mcommit b213cd503fea164f16a06b8d5263a77089469e67[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Dec 15 09:17:49 2015 -0500

    added drug_use input file

[33mcommit 85410809d5f6904f3f35ea51c083e93c18cbd060[m
Merge: 62cf6f9 a399635
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Dec 11 17:08:48 2015 -0500

    Merge pull request #38 from FRED-PITT/drug-use-model
    
    Drug use model

[33mcommit a399635e59f8472ed7080fd30a80dfea26dfb2eb[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Dec 11 17:05:34 2015 -0500

    Update Markov_Epidemic.cc

[33mcommit 6df88af2f00c2d40e16510712d057360a4428b0c[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Dec 11 14:59:51 2015 -0500

    Markov drug use model

[33mcommit 383f8c477f90a171c7ef2fc71ca39dd4f132d91c[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Dec 9 14:25:05 2015 -0500

    markov transition model

[33mcommit 74212f8910aa734c716b5cff2aad5547bc428f41[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Dec 9 12:24:31 2015 -0500

    initial check in

[33mcommit 8014ac02f8b717923dfd6399037196e6e584ac5b[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Dec 7 15:06:28 2015 -0500

    Initial check in

[33mcommit 62cf6f93a6b32c3440ac8ff06d1d826e309eed11[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Dec 4 09:23:04 2015 -0500

    update comments

[33mcommit fce26dd95092075afceec369096bd2fe1266e06d[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Dec 4 09:20:27 2015 -0500

    Initial check-in of alternative influenza disease model

[33mcommit 3434e06dc4584c5929b59e1f14bcd79b6da01e69[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Dec 4 09:12:43 2015 -0500

    update comment in defaults file

[33mcommit 1559decdee639d9d2c4976fef10a4d6ecfaacfa4[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Dec 3 13:03:30 2015 -0500

    Bug fix: incorrect asymptomatic infectivity.
    
    Also:
    
    - better support for time-varying infectivity and symptoms.
    - made binary symptoms and infectivity the default
    - recalibrate
    - new RT

[33mcommit bd76766abb9cba07d4a898ab5501e5770726101c[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Nov 19 12:46:11 2015 -0500

    Merged guido-vector-branch into master.
    
    No change in rt.

[33mcommit 39b747961e011e34f0287f487997ac9ef7b42c39[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Nov 19 09:54:07 2015 -0500

    Optional health charts
    
    Optional parameter:
    
    enable_health_charts = 0
    
    No change to rt.

[33mcommit a63a7f312a5f16a974fe36e82df260e78abd91ce[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Nov 18 11:37:45 2015 -0500

    fixed bug in disease fatality code.
    
    No change to rt.

[33mcommit 981cd396d1bf73b5748263a0c3753bfd0c861c94[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Nov 11 20:42:05 2015 -0500

    Add <disease>_make_all_susceptible parameter
    
    If true, make everyone susceptible by default. Otherwise everyone is
    unsusceptible by default.
    
    Default is true.
    
    No change to rt.

[33mcommit facb59ee309234cc32070a10098bd8cf5d6ebd50[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Nov 11 20:03:16 2015 -0500

    renamed variable Health::infected_in_mixing_group
    
    No change in rt.

[33mcommit 97c60b6a42c13b7a8b517adc0905a77c7032c29a[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Nov 11 18:07:28 2015 -0500

    Improved status messages in Health
    
    Messages ca create a health chart for each individual, showing changes
    in health status.
    
    No change to rt.

[33mcommit fed3f96998a7d9667a42d7f2d55c89f3cd15996a[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Nov 11 16:17:45 2015 -0500

    Recalibrated after last bug fix
    
    New RT results with recailbration

[33mcommit ac091a9a8dd9b513967429403e399235abdf6538[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Nov 11 15:46:57 2015 -0500

    Bug fix: symptoms now resolve

[33mcommit c08e6f205cc67cf2f8cccb501d471ed9c1c03e6b[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Nov 11 15:17:09 2015 -0500

    fixed bug that prevented people from become symptomatic

[33mcommit 61cdb5ece57717c70fddfcf561eeb34c3e4755a3[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Nov 10 20:29:50 2015 -0500

    Fixed indentation

[33mcommit 7db872ebc5940760907649c12d90b9e99865b27e[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Nov 10 12:57:16 2015 -0500

    Recalibration and new rt
    
    New rt results

[33mcommit eef9cf89501e8e7a8172d848213f9a4bbe1b6e2a[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Nov 10 09:16:22 2015 -0500

    Allow 3 digits of precision

[33mcommit 40fd5d78abc200d1df45f131f323eb7cc807fb7f[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Nov 9 18:00:56 2015 -0500

    Add fred_get_age_plot

[33mcommit 737ae214ee227cab8788e34213643b760ed25847[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Nov 9 16:16:01 2015 -0500

    Adopt age_bias transmission in neighborhoods as default
    
    New rt results.

[33mcommit 267afbc587a5307e3aac51110d011b4882668806[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Nov 9 11:51:32 2015 -0500

    Use lognormal incubation for influenza by default
    
    New rt results.

[33mcommit 3c57aa4292fb823f76efbf83a6e1997d081f3689[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sun Nov 8 16:38:29 2015 -0500

    Code cleanup and improved comments

[33mcommit c7f691b38434e6f1e0ae3f56571969c6d99ad377[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sun Nov 8 16:29:52 2015 -0500

    Make all disease-specific age maps optional.
    
    All age maps are assumes empty unless given in the params file, with
    one exception:
    
    if age_specific_prob_infection_immunity is not specified in the params
    files, we assume that INFECTION PRODUCES LIFE-LONG IMMUNITY.

[33mcommit 8f31804526ded9f69d2ff9f6d8db73030483f2a9[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Nov 7 16:59:32 2015 -0500

    Restored use of asymptomatic_infectivity multiplier
    
    The asymp_infectivity parameter was no being used.
    Restore to old usage.
    
    RTs updated.

[33mcommit dc85d73ed16a712eeefaece20d16a5a342439453[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Nov 6 10:07:43 2015 -0500

    Optional upper bounds for natural history distributions
    
    New optional parameters
    xx_incubation_period_upper_bound = 0
    
    xx_symptoms_duration_upper_bound = 0
    
    xx_latent_period_upper_bound = 0
    
    xx_infectious_duration_upper_bound = 0
    
    These only apply to lognormal distributions.
    
    If the upper bound is > 0 and the draw from the lognormal exceeds the
    upper bound, and random number between 0 and upper bound is selected.
    
    No change to rt.

[33mcommit 352f2bc3aa71d867acaf8d9202ef3fdf613941c8[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Nov 5 17:07:48 2015 -0500

    Fixed bugs in HIV_Natural_History

[33mcommit e65db816282fa6bd51de98e58dc3fa0d0d031bd8[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Nov 5 15:33:16 2015 -0500

    Reduce log messages.
    
    No change in rt

[33mcommit 08b92313307b0afbadf9a1150a94e68b320714ca[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Nov 5 10:16:52 2015 -0500

    Improved error checking in Natural_History.cc

[33mcommit 38b3b0c526613566266f63bf30191a02bc4d5968[m
Merge: 79515e9 e231331
Author: David Galloway <gallowdd@hotmail.com>
Date:   Wed Nov 4 10:34:16 2015 -0500

    Merge pull request #36 from FRED-PITT/Network_Update
    
    Network update

[33mcommit e231331460ea48ffc037408cf757311750170f9b[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Wed Nov 4 10:33:34 2015 -0500

    Update for Mixing_Group
    
    Made enroll and unenroll virtual in Mixing_Group to make sure that the
    inherited class overrides get called

[33mcommit eda95a1cf9f84e8d2e942ab0079d9a99e6ddf1a9[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue Nov 3 13:54:56 2015 -0500

    More Refactoring for mixing group

[33mcommit 79515e948c956ba6258501b9aa99bafa1ec1ca0a[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Oct 30 10:00:10 2015 -0400

    Bug fixed to sexual transmission
    
    No change to rt

[33mcommit 4a1db309b48f7e20e1cd46a872d10e8ccb761490[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Thu Oct 29 14:47:28 2015 -0400

    Updates For Networks
    
    1) Created Sexual_Transmission_Network as a subclass of Network
    2) Place subtypes are now just characters (since that is how we used
    them anyhow) and they are defined in the base class where they will be
    used (e.g. Place.h or Network.h)
    3) Moved some methods down to Mixing_Group
    4) Made some member variables in Mixing_Group private so that the
    descendants need to use the setters/getters

[33mcommit 29e0937717e810f2b2ee63043839544c3ea34e38[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Oct 27 10:00:02 2015 -0400

    Make all age maps optional
    
    If an age_map is not found, it will be set to empty.
    
    This relieves the burden on new disease definitions to provide
    parameters for all age maps.
    
    No change in rt.

[33mcommit 190a76b203095d0aefe3512c49479de5e5b0dc60[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Oct 26 16:48:52 2015 -0400

    Added optional parameters to ramp up symptoms and/or infectivity
    
    No change to rt.

[33mcommit 353295727062e62fe921064d878d7475d88be9e5[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Oct 26 08:37:36 2015 -0400

    Updated gitignore

[33mcommit 19ba84ed45bbd53f2767d08f8e1d57e010eb78cf[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sun Oct 25 09:34:39 2015 -0400

    Code cleanup
    
    No change in rt.

[33mcommit d4e7e8555b43489068d5db8bdf2b87bfab72dce4[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sun Oct 25 09:17:46 2015 -0400

    Code refactoring -- removed functors
    
    No change in rt.

[33mcommit 5ee9e593fae2da79afa09e3b29ca3b6510bfada1[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sun Oct 25 09:02:55 2015 -0400

    Minor code refactoring
    
    No change to rt.

[33mcommit f289879ad9769771e4c80410c0f72bde49c56868[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Oct 24 22:29:44 2015 -0400

    Removed obsolete code
    
    No change in rt.

[33mcommit a4b3b18db112cf80ef157cff4a43a9f4b589853d[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Oct 24 22:18:07 2015 -0400

    Removed obsolete code
    
    No change in rt.

[33mcommit 120b1784ea403598b62375d007edb1f075e94764[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Oct 24 22:10:51 2015 -0400

    Removed obsolete code
    
    No change in rt.

[33mcommit 9976c75c5a770b6f7c527ac3f1330d751483f5eb[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Oct 24 14:40:08 2015 -0400

    comments added to defaults
    
    No change to rt.

[33mcommit bf892db8684bc235b1f1eeea370cecf246c1b905[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Oct 24 12:40:59 2015 -0400

    Integrated lognormal and cdf options for infectious and symptom distributions
    
    No change to rt.

[33mcommit 39742b515cdd2d8d72a87d554af6ef4d264fbdee[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Oct 23 15:58:48 2015 -0400

    Flu parameters using lognormal distribution for incubation period and infectious offsets from start of symptoms
    
    No change to rt.

[33mcommit 781c6c8e1f6be689ceffe4a19f55491fb542bd95[m
Merge: e05ecca 1652f51
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Oct 23 09:46:54 2015 -0400

    Merge pull request #35 from FRED-PITT/incubation-offset-branch
    
    Incubation offset branch

[33mcommit 1652f51b088e744632cdeb4374d53f64cfc1afe2[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Oct 23 09:43:48 2015 -0400

    Working version of incubation_offset

[33mcommit ab0818af5fbd8ab577eab77221ca9f80ac633551[m
Merge: a50b1cc e05ecca
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Oct 22 12:36:59 2015 -0400

    Merge branch 'master' into incubation-offset-branch

[33mcommit e05ecca5a3bb69a96180c9ebf12dd77bb3c08341[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue Oct 20 13:22:35 2015 -0400

    Refactoring and Code Cleanup
    
    1) Refactored Network and Place to inherit from base class Mixing_Group
    2) Removed the unnecessary setup() from Places so that they now call
    the base constructor before they do their own initialization
    3) Code clean-up
    4) added a new parameter-based workplace size model so that we can
    modify the number of size groups as well as the max of each category

[33mcommit a50b1cc32b55a666baecb64a2c2291e2ef85155e[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Oct 15 18:06:52 2015 -0400

    Merge guido-vector_branch

[33mcommit ef2e16e90bd5dfee437ef244a2fd0886dacfacfe[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Oct 14 11:20:54 2015 -0400

    Added Date::get_sim_day()

[33mcommit c7ecc907177c9d2db538bcfa59b4bc63cae74c60[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Oct 12 16:41:01 2015 -0400

    changes to transmission age bias
    
    experimental.
    
    No change to rt.

[33mcommit 8ba3a101b0f52400e247adbe3cc3614329935175[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sun Oct 11 09:04:17 2015 -0400

    Added derived disease_specific epidemic classes
    
    Example: HIV_Epidemic.
    
    Use the derived epidemic class to maintain and print disease specific
    epidemic statistics.
    
    No change in rt.

[33mcommit 033fe6641faffeb6427453923e31cfb9d8e7863a[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Oct 8 12:42:26 2015 -0400

    Support for generating random networks
    
    No change to rt.

[33mcommit 916b309f46945fc63d931287597ea0e1bd603583[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Oct 7 17:47:56 2015 -0400

    Transmission_Network added

[33mcommit 5e08f9af05da6013756945cd17e77369b93d9c3f[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Oct 7 13:28:17 2015 -0400

    Initial support for networks

[33mcommit e1ce203ebd74df5d317493072255e109951eb5c7[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Oct 1 15:17:55 2015 -0400

    Support for same_age transmission bias
    
    No change to rt.

[33mcommit 5dccd0ddc6950e5fca458d2821dfdaf768fea660[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Oct 1 13:01:23 2015 -0400

    Initial version of fred_get_ages

[33mcommit 391bf6904facf15d5918683c1cab3179856d2477[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Oct 1 12:53:48 2015 -0400

    updated to new parameter names

[33mcommit 95b146bd4ebbe9dd9d82a7dfe0b181f15568e319[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Oct 1 09:06:38 2015 -0400

    Changed "favorite_places" to "daily_activity_locations"
    
    No change to rt.

[33mcommit 936245f33e4f77a380a0b9cf056385b86852bb99[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Sep 28 12:33:16 2015 -0400

    Changes to age map parameters
    
    Age maps parameters are now:
    
    xxx_age_groups = 1 120
    xxx_values = 1.0

[33mcommit c5937caabfa6b8dc64247b28c12488e5391156af[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Sep 28 09:31:36 2015 -0400

    at_risk parameter is now optional
    
    The <disease_name>_at_risk parameters are optional, and are required
    only if
    enable_vaccination = 1

[33mcommit 0e889fe7530246a1f551887a063be71b676c7ad0[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Sep 23 10:35:35 2015 -0400

    Added disease specific end_of_run method for producing reports.

[33mcommit f6f3999be36d134af0bc82591e1ad0d681cc8937[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Sep 22 20:33:54 2015 -0400

    memory optimization for events
    
    No change in rt.

[33mcommit 4eee9887c85e5f8e2e1dd5597596b19060d4a41e[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Sep 22 16:58:04 2015 -0400

    Cleanup up Events class.
    
    Eliminated function pointer for handlers.
    No change to rt.

[33mcommit 797b0128e774e16cd96b47534cc5a9f0065faac5[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Sep 22 14:26:40 2015 -0400

    Restored case fatality handling

[33mcommit 52f597efd555aa3e78cca3389ec02e53bc527c2d[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Sep 21 14:58:18 2015 -0400

    Renamed "trans" param to "transmissibility"
    
    “trans” takes precedent for now, but is deprecated.

[33mcommit 4991605c26c66fcc76d4e53645d5f701c872c880[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Sep 21 13:59:39 2015 -0400

    Generic respiratory contact rates
    
    Place-specific contact rates and transmission probabilities are now
    generic for all respiratory diseases.
    
    No change to rt.

[33mcommit e6cf53eb6b7d895519ade3b0dfb480d1e54575e6[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Sep 16 10:59:56 2015 -0400

    updated vaccine regression test

[33mcommit 2dec30f0aad00bba7a078f4ca58a8de01f07a332[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Sep 15 11:43:43 2015 -0400

    Update tutorial params file

[33mcommit 09ad500012f6c3dbf72f3a2579f09d0f5828f2a1[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sun Sep 13 17:06:31 2015 -0400

    update tutorial

[33mcommit 74d757358ebb9ea27181ab066bc01d94e7881e0f[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Aug 28 11:25:57 2015 -0400

    debugging statement for vector transmission

[33mcommit cab606c98d2ac11ce6d79807c8cb019fda085cbf[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Aug 27 13:28:26 2015 -0400

    Fixed bug in Vaccine_Manager and added missing parameter enable_residual_immunity_by_FIPS

[33mcommit 378fb17c020b8eae2c3c1c6339cc6a8b9780f1e7[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Aug 25 14:10:57 2015 -0400

    Removed obsolete parameter template

[33mcommit 16e30625285655b37bf5e74a3670032828d6bfa4[m
Merge: 3a1432e 62fb5dc
Author: John Grefenstette <gref@pitt.edu>
Date:   Mon Aug 24 21:22:05 2015 -0400

    Merge pull request #33 from FRED-PITT/update-neighborhood_patch
    
    Store places located in and attended by residents in Neighborhood Patch

[33mcommit 62fb5dc0fa41b262d82194f447aded46c3f0877e[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Mon Aug 24 21:19:17 2015 -0400

    Store places located in and attended by residents in Neighborhood Patch

[33mcommit 3a1432e9225b7a6b10a0a52d2283dd3145f4a726[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Aug 24 16:48:27 2015 -0400

    Bug fixes in vector disease code

[33mcommit 0a5f16252dbe0c986fc92ca23b618c3ff0e6b160[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Aug 24 14:52:01 2015 -0400

    Code cleanup.

[33mcommit 4c68d489c27c8c9d87edb5dfd6c7bc6e5eb0f6b4[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Aug 24 14:40:13 2015 -0400

    Exclude params file from repository

[33mcommit e76ff32c9f7d35b4715b9f83f6128227a58e8a47[m
Merge: b88691b ec8dcee
Author: John Grefenstette <gref@pitt.edu>
Date:   Mon Aug 24 14:29:57 2015 -0400

    Merge pull request #32 from FRED-PITT/vector-model-update
    
    Vector model update. Closes #31. Also provide stubs for Sexual Transmission mode.

[33mcommit ec8dcee07c354667878260380b050ea8aa568859[m
Merge: 509dfc2 b88691b
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Aug 24 14:24:34 2015 -0400

    Merge branch 'master' into vector-model-update
    
    Conflicts:
            src/Health.cc

[33mcommit 509dfc2bfb19a613af38621a5f44335a433863f7[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Aug 24 14:21:42 2015 -0400

    Support for sexual transmission

[33mcommit b8d84bd1ecbeaff157dfa05cb56a13b4a307ed9b[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Aug 24 14:15:43 2015 -0400

    Stubs for sexual transmission

[33mcommit 6411bd9d9a2835bc268c820deef4f6db8e52280a[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Aug 24 14:07:45 2015 -0400

    Read static parameters for place classes in Place_List

[33mcommit c5678db44905d784e7ade851f30cb2653b35df95[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Aug 24 12:50:16 2015 -0400

    Eliminate reading respiratory contact params for non-respiratory diseases

[33mcommit 1dc2280d8b5ea74205569c4f2417743a417d94d9[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Aug 24 12:08:24 2015 -0400

    Working version, but no vector control

[33mcommit 5e507a38939335ddfe0f381f165ce520b726d654[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Aug 24 10:16:37 2015 -0400

    Added vector_disease_data_t to Global.h

[33mcommit d194e8e30f819b52afc6a9258db51080b2b13dfb[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Mon Aug 24 08:42:03 2015 -0400

    Added log statements for timing

[33mcommit 63fe55f861a8bb58bcd87ec2f0ac4f048890056d[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Aug 22 13:12:25 2015 -0400

    Code cleanup

[33mcommit 9d5c277dc75f5dcf94958225578e561bb760aca0[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Aug 22 11:53:56 2015 -0400

    Fixed access to places by type

[33mcommit 7dcec3095bbc490bced0d9004ebcf0302632920e[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Aug 22 11:23:39 2015 -0400

    Cleaned up access function in Place_List

[33mcommit 150588cf2f02f88c9df1529fac2b9a6fec845eeb[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Aug 21 15:31:40 2015 -0400

    Move non-place-specific data and code to Vector_Layer and removed from Place

[33mcommit b88691bb3449da11eac466d9a4fb9c15e4637e0a[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Fri Aug 21 11:13:05 2015 -0400

    Commented out inadvertent debug
    
    I don’t think we want to have this output at verbosity 0

[33mcommit 1c0c2d23b6633194dfcc2764a8d467e4426aa58f[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Fri Aug 21 11:06:06 2015 -0400

    Quick fixes for LINUX
    
    Added <climits> so Linux compilers know what INT_MAX is.
    Captured return from fscanf in Vector_Layer.cc so Linux doesn’t throw
    warning.

[33mcommit 8907cd2e4b86fb7b101da43f64a24ccb266296f9[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Aug 21 08:59:27 2015 -0400

    Added missing file Vector_Transmission.h

[33mcommit 1cc58120d096a4887d0b9329269295fddbcfae1e[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Aug 21 08:20:02 2015 -0400

    Eliminated obsolete param visit_home_neighborhood_unless_infectious
    
    This parameter is assumed to be true now, to greatly enhance FRED’s
    runtime efficiency.

[33mcommit f9251a1a6dbc3bb0fb083b9f22be0e55742e9d9b[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Aug 20 17:35:56 2015 -0400

    Added derived classes Respiratory_Transmission and Vector_Transmission
    
    Closes #30

[33mcommit adb7c2c5047f44bfff17b2a6664325d4dbbd6d7d[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Aug 20 11:46:56 2015 -0400

    Commented out FRED_CLANG_FILES by default
    
    As noted, this can also be set using an environmental variable. This
    might be preferable to make the Makefile more portable.

[33mcommit 59f474c89d376f8cfd2be7f75a86b17967d98184[m
Merge: 565f492 9b331fa
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Aug 20 11:29:37 2015 -0400

    Merge pull request #22 from FRED-PITT/event-list--test-merge
    
    Closes #13, closes #23, closes #24, closes #25, closes #26, closes #27, closes #28, closes #29

[33mcommit 9b331fa414520375c3160b361eac08c1547c86c2[m
Merge: 3d3ff4e 565f492
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Aug 19 12:19:01 2015 -0400

    Merge branch 'master' into event-list--test-merge
    
    Conflicts:
            src/Activities.cc
            src/Activities.h
            src/DefaultIntraHost.h
            src/HIV_Intrahost.cc
            src/Household.cc
            src/Place_List.cc

[33mcommit 3d3ff4eb29d66f41ed0060d79c4c5eed7643a048[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Aug 19 12:01:39 2015 -0400

    updated comments

[33mcommit a13094657c58e0bc5e91bdd3dea010919663889d[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Aug 18 15:03:05 2015 -0400

    Support for vector transmission

[33mcommit edca72e806a34dd794b32dd842021790aee2eba7[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Aug 17 17:31:44 2015 -0400

    support for Vector Transmission

[33mcommit ea0beaaf704caee3aedfeda9d189a7625f3b8acd[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Aug 17 13:20:06 2015 -0400

    deleted obsolete code

[33mcommit 43ac53f5300a3693c52776470d0a28ada90d4541[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Aug 15 12:29:46 2015 -0400

    Eliminated static arrays in Place.h

[33mcommit 97374c6f9b0498d4b70013228a6fc3cdd30728e7[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Aug 15 12:21:05 2015 -0400

    enable_viral_evolution = 0 by default

[33mcommit 1111a57746b9bfabfc6ccdfbf7f3f3bf7ff75752[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Aug 15 11:25:11 2015 -0400

    Code cleanup

[33mcommit ba00818481192ca1e6fc2749212613117098001b[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Aug 14 15:39:02 2015 -0400

    New regression tests

[33mcommit 848cd889b3d5cf770b2debee5f41e93f21581d09[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Aug 14 15:36:55 2015 -0400

    Recalibration

[33mcommit a38c01362c7244459f1f511cc82324a6befa8a18[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Aug 14 15:01:08 2015 -0400

    Shuffle indexes, not infectious and susceptible lists

[33mcommit e35ebde86a8a59fa65c26e1c0006e604eddea6dd[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Aug 14 13:12:05 2015 -0400

    removed obsolete code

[33mcommit 992e5ebd9f70af4b57c8e3d75b0dd8ec93ca425c[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Aug 14 12:54:51 2015 -0400

    Removed disease transmission methods from Place class

[33mcommit f89e4696d3612c9f7a2bc9fc48bc2620d35328a0[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Aug 14 12:37:29 2015 -0400

    Working Transmission class

[33mcommit e8ed25f6f531a494862a5751f9cc4f1f4d40892d[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Aug 14 09:14:11 2015 -0400

    Initial move disease transmission models to Transmission

[33mcommit 1d358a800657ec3b63b8a9c3506dab61d9942ece[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Aug 12 16:29:20 2015 -0400

    Code cleanup

[33mcommit 16ef1049e9f9252dfdb775739f9a00035c5f704d[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Aug 12 13:47:45 2015 -0400

    remove obsolete files

[33mcommit bbb19d23065cf1d9725d69566f3d90e8835620ec[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Aug 12 13:39:55 2015 -0400

    Call update_infection for each active infection.

[33mcommit 30ee7c5b8e87a9462b0a425cd0b66197fbfd7f30[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Aug 12 12:12:40 2015 -0400

    Reduced status messages in Activities

[33mcommit f894ed4ab4c38a97bb7d4c5d05a516906ac99cdf[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Aug 12 12:04:33 2015 -0400

    bug fixes

[33mcommit 5d421e0fbcd252b8ba9724c02f6cc7662ba78e29[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Aug 12 11:20:19 2015 -0400

    Working version including HIV model

[33mcommit 46791c0f28136fe2678a2fcb114b9a9293440a55[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Mon Aug 10 17:54:06 2015 -0400

    Infection class finished
    
    Still need edits in Natural_History

[33mcommit 93c0f3d8ac99356b2131d91e744f15e87dccec81[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Aug 10 15:25:06 2015 -0400

    Added HIV_Natural_History class and HIV_Infection class

[33mcommit dcfa0eae853122a6052261a16a88449f611cebc8[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Mon Aug 10 08:14:51 2015 -0400

    interim edits
    
    TODO:
    Put disease-specific individual data in Infection (and derived
    classes). Put update routines in Natural_History and derived classes.

[33mcommit 565f4927f67415a22454c210b34241d663ecda46[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Fri Aug 7 15:38:34 2015 -0400

    Update compilation notes

[33mcommit 37e1f21f700f01079dd9d03b427cbcd205ea76ac[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Fri Aug 7 15:31:00 2015 -0400

    Case Change
    
    Change Intrahost.h to IntraHost.h so that the compile doesn't fail on Linux systems.

[33mcommit 6857a78fa56a8f012dca51eac555467b0e18e07a[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Fri Aug 7 15:26:48 2015 -0400

    Case change
    
    Change Intrahost.h to IntraHost.h so that the compile doesn't fail on Linux systems.

[33mcommit fff31bfb43ae4cdd52d9ef13c435e056f175c72c[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Fri Aug 7 13:38:15 2015 -0400

    HAZEL Update
    
    Changes to statistics for visualizations of healthcare deficits.
    Applies to #15

[33mcommit 84d76af734f6f04985e24d5fcb24021263f195c5[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Aug 6 14:52:20 2015 -0400

    first cut at base class Infection
    
    Still need to add methods to Disease and change name of
    Binary_Natural_History to base class Natural History

[33mcommit 9613f7b59c5bb0d8c3686e9715625f67a46e8d30[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Aug 6 11:59:04 2015 -0400

    Use dynamic arrays to reduce memory in Health

[33mcommit 0ed6ea3cdaf0f2ae5d527657ee76404893e9477a[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Aug 6 08:42:41 2015 -0400

    Undo changes to Epidemic

[33mcommit 0634888c910cc175da5db4f4f625b7a21ea1c1a8[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Aug 5 21:18:01 2015 -0400

    support for variable disease models

[33mcommit a7aac8ff2093d8ead181445c425a5dda28dd5d36[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Aug 5 16:48:50 2015 -0400

    support for HIV model

[33mcommit 116d7b2c5408665b8c8ef55d395a4abb5d6e4db6[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Aug 5 16:25:37 2015 -0400

    removed obsolete comments

[33mcommit b9ce947f8c7f648faa4680d1ffd730ac8740d4bb[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Aug 5 15:43:25 2015 -0400

    Natural_History files renamed -- no change in rt

[33mcommit da7c1ca74417a99fa16afb3930f5c47f4bf2f599[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Aug 5 15:11:16 2015 -0400

    renamed Natural History files

[33mcommit 2b605971258f7b9ca96bcd3ad7ef4e43dc983f9e[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Aug 5 14:49:10 2015 -0400

    new regression tests for recalibration

[33mcommit 2c508129193a60199ed2b1d7c34868e251b50a54[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Aug 5 14:48:02 2015 -0400

    recalibration

[33mcommit 8c9e78a2fbdcdc314cd8455dc2ca3c313d8b2af1[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Aug 5 11:49:38 2015 -0400

    Updated regression tests

[33mcommit bdaeee24d689281de947f1e91336fec22b5336df[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Aug 5 11:42:04 2015 -0400

    removed obsolete code

[33mcommit 9ba36a383c208a54effe86cb802affee31ef0b9e[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Aug 5 11:21:39 2015 -0400

    Recalibrate with contacts normalized

[33mcommit ec5da85d47a47a9ddb5a02318dd43d7c18ca328b[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Aug 4 20:55:22 2015 -0400

    epidemic timing statements

[33mcommit 4f26633adbee7419945381430a3c86ddf97634d1[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Aug 4 20:15:52 2015 -0400

    Cleanup Utils.

[33mcommit 715cf526cdb6ebbd1bed41808be41974aeaafad9[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Aug 4 17:23:54 2015 -0400

    bug fixes in high resolution timers

[33mcommit 4bb885f184bb6d12965e72bb6696fdedac3313bd[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Aug 4 16:54:09 2015 -0400

    high-resolution timers

[33mcommit a8d94fa276ecc63ea0b340f14efd8c9e9fa976b4[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Aug 4 14:44:59 2015 -0400

    working version

[33mcommit 11d0e3e1eb633bb27b52dde4384c154066a0b8df[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Aug 4 07:36:06 2015 -0400

    debugging

[33mcommit dc7b53a8d5933df589d5a2f55a6592126eb6b014[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Aug 4 07:22:58 2015 -0400

    debugging

[33mcommit 3ee3398628d01112a02542989bdf4c9b119c2494[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Aug 3 16:36:33 2015 -0400

    Compiles but no transmission

[33mcommit 2233722adeaf1c19f6368bcf3132af0ab0d64846[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Mon Aug 3 08:34:55 2015 -0400

    Removed obsolete code

[33mcommit 0de598518337129b274324fa37f3223c73511d88[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sun Aug 2 16:42:28 2015 -0400

    need to cleanup obsolete code.

[33mcommit 8e32fdff43c10191d8d39f12cab08a5ea831ddc3[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Aug 1 18:16:29 2015 -0400

    Updates to Infection class

[33mcommit 2bb833e1625cf321790d1ba38dcdfbeef2ad0fc3[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jul 31 08:45:43 2015 -0400

    check for susceptibles present in place of infection

[33mcommit 94db0304b6d39497f1a6967f454a1582796472ed[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Jul 30 22:56:37 2015 -0400

    bug fixes

[33mcommit b31bd71bc4b19b7da3dea7338a8c320081f2f8c3[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Jul 30 22:38:11 2015 -0400

    Activities updated

[33mcommit 4c8c649d1f1ef824de737cabf1af760d8795a7b5[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Jul 30 19:53:24 2015 -0400

    Redo favorite places.

[33mcommit e8510f969b5f00fefce646830ce4b35e4d3b8373[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Jul 30 09:38:52 2015 -0400

    Added Person_Place_Link.cc
    
    Cleaned up on Person_Place_Link class

[33mcommit 0cabdff0b85631700b203f1877e97dd8bb0ee5fc[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jul 29 17:50:22 2015 -0400

    Memory reductions in Infection class

[33mcommit 4a55179f910e4629f4b79c6a1e82bb5a11c08c51[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jul 29 15:13:05 2015 -0400

    interaction between epidemic and update_health

[33mcommit 1ef379583349e4fddc5eac2258276f64e2eeaf97[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jul 29 09:58:56 2015 -0400

    runs but only homes and schools infected

[33mcommit b07c11a5b2edc6dfc28d1a1ab7bfeebb005364d1[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Jul 28 09:35:13 2015 -0400

    Added Person_Place_Link.h

[33mcommit aff3eef4010be2d4b02495b5dde9a36e62bf9f15[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jul 27 17:14:53 2015 -0400

    updated active places.
    
    Compiles but aborts.

[33mcommit 2fa7ec22cc8d1a3d1c2d12c0eb0f99d96ef5b89d[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sun Jul 26 21:32:22 2015 -0400

    old code protected Global::Test

[33mcommit 043760ef90b22b7676e521c3bcacd3f0da50e8d4[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sun Jul 26 17:24:21 2015 -0400

    households and schools in sync with old code.

[33mcommit 2830167192a4bbf4653da15eefd871436cb596ae[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Jul 25 20:49:51 2015 -0400

    active sets added

[33mcommit dc1c522d14eded027a8b4ef7c5d4a8d49ca7f25e[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Jul 25 10:00:55 2015 -0400

    remove ignored files

[33mcommit 73b98398a6fb582d13c6f1de5b6151dc960ced9c[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Jul 25 09:55:54 2015 -0400

    Updated etc/gitignore
    
    cp FRED/etc/gitignore FRED/.gitignore

[33mcommit 21d12f0aff51578d52dd1408b12662ed2ce25e9d[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Jul 25 09:39:34 2015 -0400

    active_households added.

[33mcommit f07a9ee9692f78b5e6c6e17cefdf270fa71c9825[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Fri Jul 24 15:20:00 2015 -0400

    HAZEL Update
    
    Changes to Primary Care assignment. Applies to #15

[33mcommit dcf89602e1bad1b57bea6517a9f8e69b7b88a169[m
Author: mgk8 <kraulandm@gmail.com>
Date:   Wed Jul 22 14:18:49 2015 -0400

    Hazel specific map changes
    
    Added Hazel specific variables and Hazel mapping code

[33mcommit b653aa2809981ba01ead0d830b69aa0016467275[m
Author: mgk8 <kraulandm@gmail.com>
Date:   Wed Jul 22 14:16:37 2015 -0400

    Changes to map Hazel specific variables
    
    Added Hazel variables and code for Hazel specific mapping

[33mcommit 48a1822ff0e742949bc8be0ef666bcc8e53a7be8[m
Merge: d25d358 d0b52aa
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jul 22 13:06:31 2015 -0400

    Merge branch 'Refactor-school_closure'

[33mcommit d0b52aa667ab5ef66cd1da03b23e9c9d0ae40c39[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jul 22 13:05:45 2015 -0400

    Closes #19 - Refactor schools
    
    Eliminated extra container pointer in Place class.
    
    No change to rt.

[33mcommit 58c146d191abc247a94aad9e4cfa578eef25750c[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jul 22 11:45:36 2015 -0400

    Repaired school closure
    
    Fixed up parameter names and fixed bugs in school closure policies.
    School infection counts were not being set in Population::report()
    
    No change to rt

[33mcommit d25d358f533bb50eb99fab3b72bdfd6d9b161626[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jul 22 10:14:08 2015 -0400

    Removed debugging statements

[33mcommit 2e21c77beafd212586f08170e473c0f4cc86597c[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jul 22 10:08:36 2015 -0400

    Added school pointer in Classroom

[33mcommit f6ffb377080498810a6eaa6b0e26e00b05eebfa9[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Jul 21 17:02:10 2015 -0400

    Code cleanup -- Closes #20: Refactor Evolution code
    
    No change to regression tests!

[33mcommit 201daae506698c8d7241ed654a481a696d601624[m
Merge: a98eb51 a7ee048
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Jul 21 16:49:24 2015 -0400

    Merge branch 'Refactor-Evolution-code'

[33mcommit a7ee04873386aa8cff96982bc7d44fe11373fc0a[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Jul 21 16:48:25 2015 -0400

    Code cleanup
    
    Passes rt_all

[33mcommit cddfd445001870ada9c2a464568274edfe659cdc[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jul 21 12:26:53 2015 -0400

    code cleanup - no change to rt

[33mcommit 6cac22b955e6e21e0304b7b27798cd4a8233c6bf[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jul 21 10:26:01 2015 -0400

    Applied standard indentation (emacs C++ mode)
    
    No change to rt

[33mcommit fc3aedf57b5ceb94074188b0519b6b9c0abea3b7[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jul 21 10:15:41 2015 -0400

    Code cleanup -- removed Transmission.h
    
    rt_all passes.

[33mcommit 356775f5126cb15ea5d66317474958ed5bbb8717[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jul 21 08:49:14 2015 -0400

    rt_all passes

[33mcommit e877a5a4e43516f853bb7a93009a101e5561e3a0[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jul 21 08:28:14 2015 -0400

    it's pass, except AV.

[33mcommit a8a2ea706afcaf42036ba31c760fa19b1dc7fad5[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Mon Jul 20 23:09:50 2015 -0400

    compiles -- rt reports bad RR

[33mcommit 51e36d6b6a244e0a2aef6fac340be11944ec325d[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jul 20 17:20:13 2015 -0400

    Non-working
    
    Restore define get_trajectory to IntraHost.h

[33mcommit 1e44fb6fc94972cca222447bf5f4f550f8dd9d19[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jul 20 11:48:22 2015 -0400

    compiles but bombs on rt

[33mcommit af551d934583abc91a2b1ec6073b65f6c00ac3ad[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Mon Jul 20 09:01:54 2015 -0400

    Non-working snapshot

[33mcommit a9013e09a95788cac98e55c4236dd8494e54db85[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sun Jul 19 19:32:20 2015 -0400

    Non-working version.

[33mcommit 933090e6c58ce2fa199f5f7cd2c981f9adf60f37[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sun Jul 19 14:21:53 2015 -0400

    Non-working snapsho
    
    In Place()
    - decide if infection happens. this may depend on viral load of
    infector.
    - If it does have call infectee->become_exposed (disease_id, infector,
    day, place)
    - In become exposed(), create a new Infection object and record
    infector, day and place.

[33mcommit a7e87a0c49b025550a602e403d22e7323533c1bc[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sun Jul 19 13:43:30 2015 -0400

    deleted FergEvolution and supporting code
    
    No change to rt

[33mcommit 5785f472ed80cd10f49e2916fe108f98f9d7cc9b[m
Merge: d16fd45 32c0118
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Jul 18 09:01:26 2015 -0400

    Closes #17 - Merge branch 'Disease_List'
    
    Conflicts:
            src/Activities.cc
            src/Activities.h

[33mcommit 32c01181abc38c959c7a1ad804b0395d7df5e8e7[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jul 17 21:35:21 2015 -0400

    Header file cleaned up.

[33mcommit ca06cfe02ae14c92de40a37f37c954e83cbf033e[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jul 17 21:26:15 2015 -0400

    Cleaned up header files

[33mcommit ac2fe6f7423f9d462dc13096298238649ddf08d4[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jul 17 21:18:07 2015 -0400

    Cleaned up header files.

[33mcommit 4d89684103e61b8af604058b54287345d7fa1848[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jul 17 21:08:02 2015 -0400

    defined prepare_diseases() in Disease_List

[33mcommit bffe3105ed95f251f6f539c040bd79ffb21fbce7[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jul 17 20:47:53 2015 -0400

    Renamed Global::Dis to Global::Diseases
    
    Global::Diseases is now the object through which individual disease
    objects are accessed.

[33mcommit eb78045b725d41b0bac2ca3ea5ba24a12493961d[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jul 17 20:36:06 2015 -0400

    eliminated Global::Diseases parameter

[33mcommit 35ddf0288dcb9e33e1436d2f9ec6b1ad8a796012[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jul 17 20:08:13 2015 -0400

    Setting disease names in Disease_List

[33mcommit 09150bd947e0fddf0130e998d3c5284712e7f3c7[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jul 17 16:30:53 2015 -0400

    Reduce occurrence of Global::Diseases

[33mcommit 1980cf34e6f5da32ff20cb3da35a02d1d14d1b8a[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jul 17 16:15:37 2015 -0400

    Removed get_disease from Population.

[33mcommit d16fd45216b2856bc0cd21861bd55f0d692ec129[m
Merge: 3639f54 450a821
Author: David Galloway <gallowdd@hotmail.com>
Date:   Fri Jul 17 15:56:32 2015 -0400

    Merge pull request #16 from FRED-PITT/HAZEL_Update
    
    HAZEL Updates

[33mcommit 450a821b221c3b7e85eb7d67f35a1cda52e2b2d2[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Fri Jul 17 15:54:36 2015 -0400

    HAZEL Updates
    
    Changes associated with Primary Care Assignment

[33mcommit 0869d4978c1c85cfedb824cd7976d3a12cddf17b[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jul 17 15:50:08 2015 -0400

    Removed population pointer from Disease.

[33mcommit 5042083bf9f0ae58c5c5e0f7f8bc3001986d09e0[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jul 17 15:36:42 2015 -0400

    moved setup from Population

[33mcommit 9ad29f2453898397aa25375a801b07a6d2b16c7f[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jul 17 15:29:28 2015 -0400

    get_parameters moved from Population

[33mcommit 2085545a2355a813a2cc5e68b24d6a669ca43d8b[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jul 17 15:24:48 2015 -0400

    Initial version Disease_List

[33mcommit 3639f54a17f328345910c8eae029b0728ecdde9d[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Jul 16 13:01:51 2015 -0400

    Fixes #12 - Remove population pointer from Place class
    
    Removed unnecessary pointer.
    No changed to RT.

[33mcommit a98eb5197b6773fc4dfbc200e3a58f568757d9f2[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Jul 15 15:01:24 2015 -0400

    Eliminated floating point non-detrminism!
    
    FRED now produce identical results under clang and gcc 5.
    
    The solution was to use <cmath> instead of <math.h>, making cos()
    consistent across compilers.
    
    Updated rt and clanged Geo_Utils class name to Geo

[33mcommit b1c37d6d2cea5d1b1b22009af686bc10e8905094[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jul 14 23:08:58 2015 -0400

    Remove obsolete code and update rt

[33mcommit fd50c723f6b1b02873d5b72a6701129cee36ed7a[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jul 14 22:30:12 2015 -0400

    fixed type in Makefile comment

[33mcommit d50736f76b6ce7f680039bbea3d154de25a61de3[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jul 14 22:24:34 2015 -0400

    Updated rt

[33mcommit 2f7fb294a4ae55d084d79555bef8c4df295ef295[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jul 14 22:12:00 2015 -0400

    Fixes to reduce floating point non-determinism
    
    clang and gcc-5 still produce different results.
    Problems with std:sort is suspected.

[33mcommit db084fc73385e61977a11170eeb17ae9a86d7af5[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Jul 14 11:56:18 2015 -0400

    Updated to gcc 4.9.1

[33mcommit ddb027aae521670ab36d12c38ae32da30bf2406d[m
Merge: 3ffeb97 a3b746e
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Jul 14 11:42:56 2015 -0400

    Merge branch 'standard-rngs'

[33mcommit a3b746e15a147cf52807fcbf2a936923740040c0[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jul 13 14:58:19 2015 -0400

    Makefile cleanup

[33mcommit ed689ac3979a7362b92bd3d7d7d35ee0fcd86b4a[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jul 13 14:47:06 2015 -0400

    Code cleanup

[33mcommit 15344543bdc6e20e6e556f969ad08e2fb06c1223[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jul 13 14:33:12 2015 -0400

    Makefile cleanup for linux

[33mcommit 08d5faaa4fa4a7690391954e1989455d1b1fa80e[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jul 13 12:15:56 2015 -0400

    Updated seeding and regression tests

[33mcommit 6724c26c71507a6f9612a0b09820fc590c0d9108[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jul 13 11:59:35 2015 -0400

    Replace MACROs with explicit functions

[33mcommit 5e1098d0946816a54e285d063f80be2ba8d7c1b8[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jul 13 11:32:39 2015 -0400

    renamed Random classes

[33mcommit 140e7d9cfa9095ae83091df42ec00f7a0363fa00[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jul 13 11:24:00 2015 -0400

    code refactoring

[33mcommit b1ba52d357f93211589623f0cea0c537053428a6[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jul 13 10:37:02 2015 -0400

    code refactoring

[33mcommit 72071197f342300218cf7b0c824b58e3e1fc33fc[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jul 13 10:36:21 2015 -0400

    code refactoring

[33mcommit a3fcfe7fed42771c0c1b31b5dd0779632aab6232[m
Merge: 1d3d482 cf33eeb
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jul 13 10:17:12 2015 -0400

    code refactoring

[33mcommit 1d3d482d70cb14b4908e2335dc2b66e5719a7938[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jul 13 09:53:09 2015 -0400

    Merge remote-tracking branch 'origin/standard-rngs' into standard-rngs
    
    Conflicts:
            src/Fred.cc
            src/RNG.cc
            src/RNG.h
            src/Random.h

[33mcommit 3ffeb97bf84efea78e543f2b8130b95a1f2c0fd2[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Wed Jul 8 14:30:02 2015 -0400

    Bug fix for fred_calibrate
    
    Changed >& to &> so it would not get “fd error” on Linux

[33mcommit cf33eeb4cd48c728780a89a0c978f39d0b2eb09f[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jul 8 07:33:52 2015 -0400

    Regression test updated

[33mcommit 8c1ea15949659f90457e265e465b61be6e6280b8[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jul 8 07:24:55 2015 -0400

    Added irandom()

[33mcommit dc13f5ae69e97bee461e2800b7222fa60a06ab8e[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jul 7 19:35:37 2015 -0400

    Use builtin distributions

[33mcommit bde90848c9296132af1ee56a8d88d52073c8e30d[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jul 7 17:13:43 2015 -0400

    Capitalize all MACROs and remove PRNG directory

[33mcommit ff7c17de7b9580424e3f8f271b8b020ad33e8a17[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Jul 7 15:59:04 2015 -0400

    Added RNG source files

[33mcommit d357fe12fb62fbc062515ae39da7803624d85414[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Jul 7 14:54:57 2015 -0400

    Makefile cleaned up

[33mcommit 9382f1b0c41ecfa4e61692b9087b0ee1065e1bfd[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Jul 7 14:13:00 2015 -0400

    Working version

[33mcommit de88af42a9e211550123dd467da12d951d5a93b8[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Jul 2 15:37:10 2015 -0400

    switch to standard mt19937 engine

[33mcommit 69a0dfd10d2f0bd3af3684a067044585801f3fa1[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Jul 1 11:50:17 2015 -0400

    initial version of fred_get_subpop
    
    Creates a subpopulation from a list of fips codes.

[33mcommit 2b0308675712258531f762959b118f0dc261d8fd[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue Jun 30 10:16:07 2015 -0400

    HAZEL Updates
    
    Multiple Changes for HAZEL simulation

[33mcommit d53e5529cfc15bc925bfcc8ddfcda4db8d22e464[m
Merge: d48d83b 37bb049
Author: David Galloway <gallowdd@hotmail.com>
Date:   Thu May 21 10:29:37 2015 -0400

    Merge pull request #10 from FRED-PITT/HAZEL_MAY_21_2015
    
    HAZEL changes and merge code

[33mcommit 37bb0493e5de00c21edf2203a040a3fec851736a[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Thu May 21 10:28:56 2015 -0400

    HAZEL changes and merge code
    
    Some HAZEL changes and merge code with new County code

[33mcommit d48d83b802ce8b4c996342fa66f38055475ecbc2[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed May 20 16:09:23 2015 -0400

    Updates to fred_make_movie

[33mcommit 0184a2a794e793f95b445281ee35ee02472268ac[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed May 20 14:22:34 2015 -0400

    Put gitignore into FRED/etc

[33mcommit d04e62afd388a30cc3cddfc47d064f74798d10bf[m
Merge: 0ae6a57 fda2e34
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed May 20 14:17:40 2015 -0400

    Merge branch 'Demographics-dynamics-May-2015'
    
    Conflicts:
            src/Vector_Layer.cc

[33mcommit fda2e34e1449a5dd5ebc184159b27a32b8e117e2[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed May 20 13:20:14 2015 -0400

    death_rate_* renamed mortality_rate_*

[33mcommit 0fa74c0f366eda3f53e862a0ee6c49e80dbae980[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed May 13 16:12:43 2015 -0400

    updates to County class

[33mcommit 95556126f214eb567f4cf6f6102c423526d065f8[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue May 12 20:55:16 2015 -0400

    Update rates semi-annually

[33mcommit 295632c42df4254596bf89d1dc9d08ce8d884487[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue May 12 20:19:22 2015 -0400

    Add County class

[33mcommit 790fe18ad3269b6f1a40d84ed5ddc3d5763d6a43[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue May 12 17:25:44 2015 -0400

    Working version

[33mcommit 0ae6a570a1fc1dcb128b6e021fa02107e4879fc1[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue May 12 10:36:01 2015 -0400

    Changed rt
    
    Update to rt

[33mcommit 3834cb6e8082b50106a44f4c192c9cb29d571796[m
Merge: d5e92a3 fb72eb7
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue May 12 09:23:54 2015 -0400

    Merge pull request #9 from FRED-PITT/HAZEL_Update
    
    Updates for HAZEL

[33mcommit fb72eb7504022a300823137e9f1981da2ffae4c2[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue May 12 09:23:24 2015 -0400

    Updates for HAZEL
    
    Updates for HAZEL as well as some formatting. Added additional fields
    to Person::to_string()

[33mcommit 7f79397de46165b92583665323da91f40ebbed22[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue May 12 09:10:59 2015 -0400

    Initial refactor for new County class

[33mcommit acd6fda55454a3e10967bfaaa9598349b75f03d6[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed May 6 15:25:39 2015 -0400

    Update on July 1, to match census data.

[33mcommit f5ecf54f68e3121a5be22aa201e3fc9d4df682a8[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed May 6 14:08:19 2015 -0400

    Population size control via adaptive death rate

[33mcommit d5e92a307f68ed7715572da36530d62e7de0cdef[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue May 5 13:30:47 2015 -0400

    Add so sprintf is available

[33mcommit cd86b4f48bc9e14080eb8fac442b460627866d03[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue May 5 13:23:51 2015 -0400

    Formatting and bug fixes
    
    Added some bug fixes to remove warnings when compiling on linux

[33mcommit 79dbfa16871d3fa21c1fde32ab2f473d6f86f01b[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon May 4 16:30:03 2015 -0400

    Compensatory death rates

[33mcommit 21aeeee997a6bb56580419e2075eb826b2975416[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Apr 30 13:19:52 2015 -0400

    restore files mistakenly deleted

[33mcommit 924e033a214b7da3794995e549726f1701ae7539[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Apr 30 12:40:37 2015 -0400

    update demographics dynamics

[33mcommit 189d3f72fa9d2d5a63674ea40a12f8282946a8e5[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Apr 28 18:33:08 2015 -0400

    Fixed bug in county population counts

[33mcommit 71a7817eef18299913edb80b0a911d352c68990a[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Apr 28 14:46:11 2015 -0400

    Support for seasonal reduction in transmissibility

[33mcommit 8da0daaa061ebd55dda4ef4ff79ee953fe8abb96[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Apr 23 16:57:51 2015 -0400

    Support for calendar-based seasonality

[33mcommit fc404fbed41a1ada7389ed6c774dd04e60c921ea[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Apr 23 14:33:31 2015 -0400

    Added mortality_rate_multiplier param

[33mcommit 8e54a008efdb6f10e8d565f057bdceb77aa7b3a7[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Apr 23 14:32:07 2015 -0400

    Default params added

[33mcommit 9a95f44c94dcc69a3268ea64d74d7d481128327f[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Apr 23 14:20:52 2015 -0400

    Added --finish and --point_size options to fred_make_movie

[33mcommit 23775e03e08f54833fd97277a61927ba9c8ddf60[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Thu Apr 23 14:14:09 2015 -0400

    Update to Susc Mod
    
    Altered susceptibility modifier to be linear up to the 90th percentile
    of HH Income rather than to the max. Any HH higher than the 90% will be
    modified at the same rate.

[33mcommit af085aaac4cf04dce904643fb002b0d7d41cffff[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Apr 23 14:13:11 2015 -0400

    Added --lw linewidth option to fred_plot

[33mcommit 9c5064d6cabffd6290813f91ea809b0672221dca[m
Merge: adb76f7 b52bbf6
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue Apr 21 12:52:38 2015 -0400

    Merge pull request #8 from FRED-PITT/ExperSuscByHHIncome
    
    Update function

[33mcommit b52bbf670306c0e464a4ff5c123e153e8a635157[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue Apr 21 12:51:31 2015 -0400

    Update function
    
    Updating function for experimental susceptibility modifier due to
    household income

[33mcommit adb76f7e4476abfbacb3aca5b408aebd2fdfbf42[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Apr 21 11:56:51 2015 -0400

    Support for birth_rate_multiplier restored

[33mcommit 5ee7512bf1f11173558bfdd5e88f64f848bafe6c[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Apr 21 11:53:35 2015 -0400

    Cleaned up rt test files

[33mcommit 941cf63a1dbc44af4f65e4efe02fc5b552cf7ff9[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Apr 21 11:39:25 2015 -0400

    Delete obsolete test files

[33mcommit 267665fce63c5c6e96b78b7aae52890cbcd6b3d1[m
Merge: e9da4cf 89e8e36
Author: David Galloway <gallowdd@hotmail.com>
Date:   Mon Apr 20 09:31:34 2015 -0400

    Merge pull request #7 from FRED-PITT/PartialChildhoodPresenteeism
    
    Update of rt

[33mcommit 89e8e36b30522c344a6ec61ed5907a40d25a4705[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Mon Apr 20 08:10:15 2015 -0400

    Update of rt
    
    Changes to Activities means rt change.
    Change to red_plot to allow for png to be created if wanted instead of
    pdf

[33mcommit e9da4cf755e28e15328299975f997daaff20e754[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Apr 17 16:48:46 2015 -0400

    bug fix to fred_make_maps

[33mcommit d33894aeb2f9a7dd4b41e263649beda5f373020f[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Apr 15 17:51:20 2015 -0400

    Support for cholopleths

[33mcommit bfc834d3c5404af030b9be1b854ccd713212b1c1[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Apr 15 17:36:40 2015 -0400

    Support for tract level movies

[33mcommit 390e8f2b0edeb3a4c6031256c35e667ac0ed7a67[m
Merge: 3dbb1c6 49bc5ec
Author: David Galloway <gallowdd@hotmail.com>
Date:   Wed Apr 15 16:45:53 2015 -0400

    Merge pull request #6 from FRED-PITT/MaryMeasles
    
    Initial Changes Immunity

[33mcommit 49bc5ec7806dc4fa13e1ee5faa0d530237692c16[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Wed Apr 15 16:45:10 2015 -0400

    Initial Changes Immunity
    
    Changes for residual immunity by FIPS code

[33mcommit 3dbb1c6bf037832c61494b63229f8b2e47e52955[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Apr 15 10:15:16 2015 -0400

    Bug fixes to fred_make_movie

[33mcommit cf66b6fbeb96fb240264e817b2c5c5ca44e035bf[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue Apr 14 11:01:14 2015 -0400

    Incremental Sick Leave
    
    Incremental Check in of childhood Presenteeism code that does not alter
    the rt

[33mcommit ad1c060f05cb69b99639a2ab11f60fed41a3d535[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Apr 13 17:27:13 2015 -0400

    Updates to fred_make_movies and fred_calibrate

[33mcommit 3bfe8035ae812001e221f4f1e70ae82b7b443a36[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Mon Apr 13 12:34:17 2015 -0400

    Back out of School Presenteeism
    
    Blackout of the School Presenteeism changes because they need to be
    more thoroughly tested.

[33mcommit cdcd1fb24f4c16b967212ce09abfff5a89fbe04d[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Apr 11 09:38:44 2015 -0400

    Report transmission by ages and upgrades to fred_plot

[33mcommit fdb4193e7b387472aabef447a7ceb8f052bcb76e[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Thu Apr 9 10:34:11 2015 -0400

    Accidental Commit Removed
    
    This file got mixed into the commit ... sorry

[33mcommit c4f9b6c98fab96daa31d0b68e4601d408c39105f[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Thu Apr 9 10:31:26 2015 -0400

    School Sick Leave
    
    Partial Changes for Childhood Sick Time

[33mcommit cd35e2abd477a802cab36c7cc7989a04b17cf0ac[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Apr 7 14:46:32 2015 -0400

    Updates to fred_make_movies

[33mcommit ad12ed1d6feed9dd304357d434e04c4de07588bc[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Apr 6 16:58:32 2015 -0400

    update fred_make_maps

[33mcommit a36ebd0a8c35698269b2d3a4f08751f31f367e04[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Apr 6 11:31:05 2015 -0400

    Updates to fred_make_maps

[33mcommit 975f0f89acb8e03e5836dd1f0061e78333fdf1a4[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Apr 2 16:46:56 2015 -0400

    Updates to fred_make_maps

[33mcommit 6923db0a81aa1d5799563601044b21d22c8e94a7[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Mar 31 12:54:52 2015 -0400

    Added api_key to get googlemap

[33mcommit a6fa3b4fe798528dab7a202216e3cdd319b57352[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Mar 30 13:26:58 2015 -0400

    Updated fred_make_map

[33mcommit ec8222037bfdf940368470ef8c3c5e24f1538723[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Mar 30 13:06:02 2015 -0400

    Updated regression tests
    
    Updates regression tests due to a change in defaults:
    use_mean_latitude = 1
    This improves FRED maps.

[33mcommit 8355937a89ba87d6055bb5ec09a9f96daf5aab2e[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Mar 30 12:41:50 2015 -0400

    Updates to fred_make_movie

[33mcommit 47d026d14566f11db7b924ffa0d025a357ebcca1[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Mar 27 21:55:18 2015 -0400

    Updates to fred_make_movie
    
    Add support for google maps.

[33mcommit 6530495403dd07f5fd45bb084359b736e7c0d795[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Mon Mar 9 10:32:35 2015 -0400

    HAZEL Additions
    
    Additional code for HAZEL including parameters to visualize the HAZEL
    outputs

[33mcommit 3d6721bf203e26d26a051c866adc16537f9325a8[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Mon Feb 23 14:19:38 2015 -0500

    HAZEL
    
    HAZEL, Code Standardization and Initial Paid Sick Days and Schools code

[33mcommit 18ddd5b03539fecbbe7550f7bf8b272fb96083e0[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Feb 12 08:33:27 2015 -0500

    Updated: fred_make_movie -h for instructions

[33mcommit 2ec6afee962d7b02f85f5a268d738c7372644763[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Feb 11 16:11:45 2015 -0500

    Updates to FRED maps and movies.

[33mcommit 6dc0107b520bbee649d98903bf01204f32bfdd2d[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Feb 11 12:45:43 2015 -0500

    Updates to FRED visualization scripts

[33mcommit 2dec6c9ca929b091ab7e01293122280d803736ff[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Feb 11 10:30:33 2015 -0500

    Removed obsolete script fred_get_county_incidence.
    
    This functionality is now available by setting
    
    report_incidence_by_county = 1

[33mcommit aa896ec6b2e702f37da3cb3615646084bda9cc20[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Feb 11 10:23:19 2015 -0500

    Allow adding new runs to an existing job.
    
    The following will add a run to an existing job:
    fred_job -k key -f -s run
    
    All LOG files are now store in the DATA/OUT directory, whether the job
    runs in parallel or in serial mode.

[33mcommit 8b9c83746272777f88d07f593ad7f8a6a48df162[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Feb 10 13:51:47 2015 -0500

    fred_cases returns number of cases in a given run.

[33mcommit 3ef1d7c50c9bd73c36c107f0046f61165601a717[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue Feb 10 09:38:41 2015 -0500

    Makefile and Code Standardization
    
    Updated Makefile to used the proper include directories for Snappy
    Updated a few files with code standardizations

[33mcommit 0e9905e6d59a8abbe366ccbce5d2d75057193bdc[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Feb 9 16:57:41 2015 -0500

    Update visualization methods

[33mcommit 42bc03c2b6831c0bbd8cd59a7ff87c5e0c1f22e2[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Feb 9 13:06:12 2015 -0500

    Visualization updates

[33mcommit 6fb2af9bf6397e400e331cd5ed21cbdbfee63ddc[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Feb 9 09:07:33 2015 -0500

    Visualization upgrades

[33mcommit aa61dd40c57054263620888529686b4fe70c36c3[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Feb 6 12:01:15 2015 -0500

    Bug fix: head of households fro group quarters. rt updated.

[33mcommit ff2c4ebc78a79f4fa57367e266bbb3cfed0fccab[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Thu Feb 5 14:56:44 2015 -0500

    Add option to not display movie file.

[33mcommit 355faa50ec2daf7f207034f99b933591d217dd7c[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Feb 5 13:57:36 2015 -0500

    Add non-dislay option

[33mcommit e10e61b27fc85b91b6446f0239f38512908f4d35[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Feb 3 20:01:05 2015 -0500

    Remove obsolete code

[33mcommit 5edfb4b3a50660aa110af6248236869beba7c8f0[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Feb 3 16:58:43 2015 -0500

    Updates to case importation.
    
    - Move processing of primary_cases files to Epidemic.cc
    - Fixed bugs in location-specific seeding.
    - Age-specific seeding takes real values for age limits.
    - Eliminated Multistrain_Timestep_Map class.
    - Updated rt.

[33mcommit e5f5c33c595f5f4bed237e4842df2bb5922fd503[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Jan 30 17:08:37 2015 -0500

    Reduced more log messages due to travel.

[33mcommit e7e5d661996a30dd7e8ff4e29a9d7154391c9d15[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Jan 30 16:49:57 2015 -0500

    Reduced log messages due to travel.

[33mcommit 778511536c7bf701116e778ead4f3593f99c444b[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Jan 30 14:38:42 2015 -0500

    fixes for new out file format

[33mcommit af26624d41f2b2e1b4e74e8e2f5e9898eea50d49[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Fri Jan 30 09:34:41 2015 -0500

    code optimizations

[33mcommit 09eebf0a769874ae491c87f2abbaa766fd8c3d44[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Jan 29 17:20:30 2015 -0500

    Use system snappy installation by default

[33mcommit 2220eaec188a59800b645bbbe8c344f1a1e17b12[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Jan 29 16:57:20 2015 -0500

    Fixed cross-compiler bug in gravity model setup.

[33mcommit f06759dfec225122fe12eeb394694770337d5085[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Jan 29 13:42:40 2015 -0500

    Handle optional tr1 library.

[33mcommit 88e441e818145ba49fa0f2e66023806dd18c9cfa[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Jan 29 11:22:16 2015 -0500

    Renamed new_master branch to master.

[33mcommit 93df1cb5d760535b8bcd344e087b6439089e6362[m
Merge: c788f77 0c41712
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Jan 29 11:15:12 2015 -0500

    Merge branch 'master' into new_master
    
    Conflicts:
            src/Activities.cc
            src/Person.cc

[33mcommit c788f7740aea4f61b501a08dd197912ae600072d[m
Merge: c1ee1c3 96c6a28
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Jan 29 10:40:32 2015 -0500

    Merge branch 'FRED-v2.11.0-Development' into FRED-event-list
    
    Conflicts:
            src/Activities.cc
            src/Person.cc

[33mcommit c1ee1c32b4ca27ce3984874349dd47d0fb464d62[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jan 28 15:11:48 2015 -0500

    Reduced log messages.

[33mcommit 96c6a28d3b136556d9b72bc5866af6259cb50eca[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Wed Jan 28 09:58:48 2015 -0500

    Hospitals Now Accept Insurance
    
    Changes to allow Hospitals to accept certain insurances

[33mcommit 97c9cb06ed4bc5b621a27ecb5013bb0fe025fb08[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jan 27 15:50:01 2015 -0500

    Initial support for travel events -- not finished.

[33mcommit 8fcf783add819a594ddcd5e63ea4fc32f0077d8b[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jan 27 14:19:37 2015 -0500

    code cleanup

[33mcommit 9b483fed4aea9c50f878eee30ee0cf265386cae7[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jan 27 14:04:41 2015 -0500

    refactored all dynamic demographics in Demographics -- no change in rt.

[33mcommit 78f279aefd8a4a84d6ad926cabe26f06e8e5550d[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jan 27 13:11:10 2015 -0500

    Handle deaths via events.

[33mcommit 924d6ff89c25d8f7389229111d432e843e028494[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Jan 27 10:50:12 2015 -0500

    Births handled via events.

[33mcommit 36202d74195e06d3ecced3b95d79b7e311ab103f[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Mon Jan 26 22:03:53 2015 -0500

    code cleanup

[33mcommit 7b70676dbd263d518367e9635a9a76f77ef54c48[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jan 26 17:15:46 2015 -0500

    Add Events class

[33mcommit 573f7f6b89e7981febb18198e010c567fc1b3bd6[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jan 26 17:09:45 2015 -0500

    debugging

[33mcommit ed61a2aa3fd306eea26abe7a2259c961c1c60410[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jan 26 16:44:09 2015 -0500

    support for event lists -- conception and births

[33mcommit 427bc7087ee608fe557ef38be1140127e6be5083[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Fri Jan 23 11:04:39 2015 -0500

    Hospital Changes and Code Standardization
    
    Added some additions for hospitalization including fixing a bug that
    was occurring if an agent had a Hospital assigned already.
    Also standardized some of the code for getting households if agent is
    hospitalized.

[33mcommit fe04cfc55c7deca8f30728b4a8826e0b7d455fdc[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Thu Jan 22 09:53:55 2015 -0500

    Hospitalization Changes
    
    Change how Hospital assignment is done so that we can take into account
    a max distance and Insurance accepted

[33mcommit e544d6e4755366e48a27bf591ff6dc2a3f0e15d0[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Jan 21 17:22:28 2015 -0500

    Code optimization for long-period runs.

[33mcommit de62e7d81250ddd8b03ccfa54eaa91f8510ddfc8[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Jan 20 17:11:07 2015 -0500

    Code optimizations. New rt results.

[33mcommit 42c35fabb5082ab718bb05193824dd9d6f201e5d[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Jan 17 17:22:28 2015 -0500

    delete obsolete file

[33mcommit 223b612af41ca758b8809519a8a857eed93e4673[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Jan 17 12:38:45 2015 -0500

    Simplified Date class. TODO: Seasonality

[33mcommit 5c091e4e17f772a15f0a9ce069ac10567b237557[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Jan 16 12:13:14 2015 -0500

    fred_plot -r option; updates to Date class

[33mcommit aa61918302ba530858617cb524dcdb153696d0b0[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Jan 15 12:31:24 2015 -0500

    Suppress printing some debugging files by default.

[33mcommit f19cf68e95971b555fab1de70ad08a0836fab184[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Jan 15 09:51:52 2015 -0500

    Bug fix on fred_job -- changed name of birth/mortality rate files

[33mcommit 295afb7b4103402b25610e9037604e8785807bab[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Jan 14 11:09:38 2015 -0500

    Support for location specific birth rates.

[33mcommit 362d8cb523e7d8308ccd4309106b0ebc9805a482[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Jan 14 10:52:12 2015 -0500

    Changes to control population growth rate.

[33mcommit ce1083fbfdad2cf0e0a8e19bbbcef7db42a15c08[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jan 12 15:54:23 2015 -0500

    Changed EPOCH_START_YEAR to 1700

[33mcommit 274e8b2685a48946aa9a116aa9888341683163ec[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jan 12 15:37:56 2015 -0500

    data file for US birth and death rates

[33mcommit c366d9abc7983a5176812294e518da46d1d3e492[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jan 12 15:09:41 2015 -0500

    deleted obsolete parameters birth_rate_multiplier and max_adjustment

[33mcommit 968d7f39594996a47666a1588bf1fb84d792664b[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jan 12 14:33:07 2015 -0500

    Support for user-selected population growth rate

[33mcommit b7ff55761fe21328695baa23ed32be77868868bf[m
Merge: 60364bb 97c187a
Author: David Galloway <gallowdd@hotmail.com>
Date:   Mon Jan 12 14:27:10 2015 -0500

    Merge pull request #5 from FRED-PITT/HAZEL
    
    Update to Health Insurance

[33mcommit 60364bb5d5479ed24af70a8ff5a3c0ee590df3b8[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Jan 12 13:17:35 2015 -0500

    Initialize target_popsize

[33mcommit 97c187a4251824add780dd8c79cf966891f12b61[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Mon Jan 12 10:49:51 2015 -0500

    Update to Health Insurance
    
    Will set everyone in a Household to the same Insurance Type
    Default insurance is UNSET

[33mcommit 75d8b840a82bab88cc1950eb48e396f36ed37d12[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Jan 9 16:56:48 2015 -0500

    updated regression tests.

[33mcommit 77fa06e3782463526baa0d02c3da3292709aea56[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Jan 9 16:44:48 2015 -0500

    Updates supporting birth and death rates.

[33mcommit e550b445c8db73269936004486b3662b2e13476f[m
Merge: 0ba2eb0 e98d9f0
Author: David Galloway <gallowdd@hotmail.com>
Date:   Thu Jan 8 12:46:12 2015 -0500

    Merge pull request #4 from FRED-PITT/HAZEL
    
    HAZEL push

[33mcommit e98d9f0f690e434cfad05ea1e8228fac8568e4b0[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Thu Jan 8 12:45:06 2015 -0500

    HAZEL push
    
    Added Health Insurance CDF

[33mcommit 0ba2eb02252e35cbf14470cbeefe9877362b0b7f[m
Merge: 952d700 921e8f1
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue Dec 23 11:56:32 2014 -0500

    Merge pull request #3 from FRED-PITT/HAZEL
    
    Code cleanup and  Insurance in hospitals

[33mcommit 921e8f126da090803d41b3a69c981a60919bc0e8[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue Dec 23 11:55:24 2014 -0500

    Code cleanup and  Insurance in hospitals

[33mcommit 952d70047651897182e44bdde1a94d76396dc86c[m
Merge: b18ab16 bef0280
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue Dec 23 08:51:42 2014 -0500

    Merge pull request #2 from FRED-PITT/HAZEL
    
    Update Copyright dates

[33mcommit bef0280237d436a5f9b49b6c7ad8c0b4689bda2f[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue Dec 23 08:49:53 2014 -0500

    Update Copyright dates

[33mcommit b18ab16361ae5aac8c19b92e849960cdd669dc98[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Mon Dec 22 15:46:23 2014 -0500

    More changes for HAZEL and Chronic Conditions

[33mcommit edfeef6d90b3c5a843559400d5e3bf33e0f69ff6[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Dec 19 14:09:15 2014 -0500

    fred_plot uses different colors for each run with -a

[33mcommit 548ee35e778b09cd2c691faef6f0606e2face208[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Dec 18 15:54:13 2014 -0500

    Fixed demographics for historical simulations.

[33mcommit 9919a82f8e874b34a96392843abe35810a115a15[m
Merge: 7867231 4bdc5a1
Author: David Galloway <gallowdd@hotmail.com>
Date:   Wed Dec 17 10:06:58 2014 -0500

    Merge pull request #1 from FRED-PITT/HAZEL
    
    HAZEL changes to FRED

[33mcommit 4bdc5a1fe5e72aeae80c0de34ab2552874c0cdbd[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Wed Dec 17 10:00:26 2014 -0500

    HAZEL changes to FRED

[33mcommit 7867231f47e83ee2857a2c4a70ab03e6239025d6[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue Dec 16 13:39:59 2014 -0500

    add default values for new age maps

[33mcommit 0fe09adad0480d185e73cba698af85c43a9f728a[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue Dec 16 13:05:19 2014 -0500

    Initial Check-in of HAZEL additions
    
    Added two Chronic Conditions and added Insurance Assignment

[33mcommit 851f6a619955ca5c43d59c3b3406d9ae5ae5a970[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Dec 12 15:11:36 2014 -0500

    Improved support for population dynamics.

[33mcommit 200907504df4b64f3476c23cf4cbfdfa3009a2e1[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Dec 9 14:23:55 2014 -0500

    Support for Colombia vector model.

[33mcommit 5877ecd48502a64b8b1860866764577fd756d140[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Dec 3 16:46:09 2014 -0500

    Fixed bug: residual immunity does not apply to newborns.

[33mcommit b4da3b1c24261f2f4ed7f388917ed88c325da414[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Dec 2 17:58:22 2014 -0500

    Fixed bug in Health::become_immune()

[33mcommit f29d3217ad8e27100158ceac93a0b61b440597e6[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Wed Nov 26 13:57:51 2014 -0500

    Updated vaccine regression test.

[33mcommit 9bacc70b2e5b6688fc27cd464a1c26016149a57e[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Nov 25 16:46:58 2014 -0500

    Support for vaccine waning

[33mcommit 0a9a11e429a14cf4c058b66475948f78b80a6444[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Nov 25 15:59:45 2014 -0500

    Initial support for vaccine immunity waning.

[33mcommit d3571ecb0fc120de8b5885db01f2ad4aad8679a5[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Nov 25 10:04:19 2014 -0500

    Updated parameter help files.

[33mcommit 79f4d067b1e0f6f6149434f5c7df1bf549980f11[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Nov 25 09:56:39 2014 -0500

    Updated fred_help command.

[33mcommit 30c43bdf831e00c7a8ed1c6a621a09bfccd9cdb7[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Nov 24 14:15:01 2014 -0500

    Updated parameter help pages.

[33mcommit a8c0d9eb7fb763dc24f95bdfd17c80c53bd26f93[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Nov 24 11:42:57 2014 -0500

    Parameter help files updated

[33mcommit 1f64270fd4c7f1b398f0f53f017ad06086a34676[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Nov 22 16:35:20 2014 -0500

    Update parameter help files.

[33mcommit fe18af9e886be44e055a4c0942c5f6b6c5f26c00[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Nov 22 16:33:14 2014 -0500

    Updated parameter help files.

[33mcommit ecfc4716d965bc7af1fecb141487c39b0f3dbb60[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Nov 22 16:19:46 2014 -0500

    update help page for quality_control parameter

[33mcommit ad908a25e14f17b9f321c86e9e1a572f6f24b455[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Nov 21 16:16:39 2014 -0500

    Added parameter documentation lookup to fred_help

[33mcommit a18ebd429de167d5affb2c8f8488affec51c802e[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Nov 21 14:06:21 2014 -0500

    Support for age maps based on real ages.

[33mcommit d8ed1800b2192008751bede17c4324f9089e275b[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Nov 20 15:15:43 2014 -0500

    Age specific probabilities of symptoms.

[33mcommit a63ed769b836bf7b5f68f01e478ade434d5364fa[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Nov 20 13:13:15 2014 -0500

    fred_sa script for local sensitivity analysis

[33mcommit ca1fe34a9fdbec2ad43b5b840b2e362472573145[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Nov 20 13:06:43 2014 -0500

    Support for plot by age per 100,000

[33mcommit 488d383bf42bd686bdcbbdac08ffba6c898d7a08[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Nov 20 12:24:09 2014 -0500

    Output incidence by age if report_age_of_infection=3

[33mcommit 7f90c1dcab9a062d787c684dee160776485fca0b[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Nov 18 20:33:28 2014 -0500

    Use clang on OS X 10.
    
    Still using deprecated -stdlib=libstdc++

[33mcommit 65559768eb2917e80d1aa204050affce8178ae21[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Nov 18 12:10:33 2014 -0500

    Makefile changes to support clang

[33mcommit f8d992f172074ae867e867f3a5cfc8e5d014ef0a[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Nov 12 11:19:56 2014 -0500

    Added -V option to fred_plot

[33mcommit 3ac1b0dfabfdcd239d583fb52cf189c17c338d6f[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Nov 4 13:57:39 2014 -0500

    Updated help message for fred_plot

[33mcommit 8830bc186c9fa680a869bf08d22306f4aac005f1[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Nov 4 13:10:47 2014 -0500

    Included disease name in infectivity profile params.

[33mcommit 9b6935938a039aa7eb4c7304a7ca86f9a28aa944[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Nov 4 12:42:19 2014 -0500

    Moves measles params to params.measles

[33mcommit 299fb725404a1058763bd2b30150b6732003f807[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Tue Nov 4 11:28:47 2014 -0500

    Add measles fixed trajectory to defaults
    
    Add measles fixed trajectory to defaults

[33mcommit bd5633b671fe6349ed83ba4e725d2d5c1fc09ca4[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Mon Oct 27 11:43:05 2014 -0400

    Comments added to defaults file.

[33mcommit 279a47ef22560f459b3b10e47d287db71c6cd4c4[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Oct 21 10:04:04 2014 -0400

    Added short option for font size in fred_plot.

[33mcommit 76fa1c7c33d47f2425119284c5ef5d26869504b3[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Oct 21 09:53:53 2014 -0400

    Add short options for most fred_plot options

[33mcommit 5f212e852bae9582e3fb7cc759324253b9ea7dd1[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Oct 21 09:15:42 2014 -0400

    Added fred_peak command, updated fred_plot

[33mcommit ed378ee7474bad17cc1b3b7f06e6a4483c383103[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Oct 17 16:09:50 2014 -0400

    utility script to generate a truncated exponential.

[33mcommit 754640db5594709507ea3d14601df80d02bd9610[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Oct 17 14:29:12 2014 -0400

    Change output file variable to Tg for generation time.

[33mcommit d64c66799734ed82205ee537005b86eb472f8802[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Oct 17 14:08:32 2014 -0400

    Updated fred_make_params and fred_make_tasks.

[33mcommit 13e742c755d62d4858e2e5c1104b9b095d1850aa[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Oct 17 12:32:45 2014 -0400

    Fixed bug in fred_plot on multiple bar charts.

[33mcommit 0b2c823719e5a358e6e15d0829f75841a0c25502[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Oct 17 12:27:49 2014 -0400

    Added a simple case isolation model.

[33mcommit f1fecc07f7bec741bf0d66e407943c303fb9c161[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Oct 17 10:56:10 2014 -0400

    Renamed shelter_in_place parameters

[33mcommit 084db215acb7c27e8556a83321ef5df2c8cc89d4[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Oct 17 09:35:51 2014 -0400

    fixed bug in fred_get_places
    
    Bug caused by recent change to outfile format.

[33mcommit 4e57c88e26e818a29452e22134287ff85aba60da[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Oct 16 17:11:42 2014 -0400

    Improvements in fred_plot

[33mcommit b74cb67757d470c07f98954a64e73718b40cb01b[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Oct 15 18:01:45 2014 -0400

    fred_plot updated
    
    Many new options provided.
    See “fred_plot —help”

[33mcommit 24a3fd9bb88d9342dc5914c376b3a147ea13ecb3[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Oct 15 10:37:12 2014 -0400

    FRED commands cleanup
    
    Improved fred_post_process and related commands to extract data from
    fred jobs. Some changes to output file format, so new regression tests
    were run.

[33mcommit 4cfcff0b80d18237d31b5de40685f3f9ed00f63c[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Oct 7 17:43:37 2014 -0400

    Added ensemble option to fred_plot

[33mcommit f148a5b2fec5812ba8b619e525f719061b715316[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Mon Oct 6 15:23:58 2014 -0400

    Updated generation time reporting.

[33mcommit fe6fc237215f204c89c0599b81f74dabb4190ca0[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Sun Oct 5 17:13:55 2014 -0400

    track_infection_events updated
    
    Added infector_exposure_date if track_infection_events > 1

[33mcommit 6d66873a5f95f2736232901f32003f1cbebb2a0d[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Sun Oct 5 16:48:43 2014 -0400

    Updated generation time reporting

[33mcommit 5a399f0235a65b7e8b06e1eb9474d80202ab0806[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Sun Oct 5 09:57:06 2014 -0400

    Updated generation time reporting

[33mcommit 70acd91f8c5c24bff8aa92a508e9b7b963cb2a2d[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Sun Oct 5 09:49:51 2014 -0400

    Updated generation time reporting

[33mcommit 19be4f05752594b7b86976442af1b98228ef758e[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Sun Oct 5 09:16:05 2014 -0400

    Shelter in place updated
    
    Changed terminology from “isolation” to “shelter”

[33mcommit 426654dab023d021cea7f3d7c199d8f855e4b663[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Sat Oct 4 09:15:08 2014 -0400

    Added etc/gitignore

[33mcommit ed4bf8a2eadfde34fb0c8a12edddcccf98ad2bf1[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Fri Oct 3 09:00:01 2014 -0400

    Changed name of params.default file to defaults
    
    The params.default file is now called “defaults”. No change to rt.

[33mcommit 1cdb57f35b202f1d7731e375ce9bfbd13222ec6f[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Oct 2 17:26:08 2014 -0400

    updated FRED visualization tools

[33mcommit 2ea9b2898bc9a604b512279f5f67d97647b06dce[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Thu Oct 2 09:31:15 2014 -0400

    Updated FRED visualization methods

[33mcommit 826c00b6fea60ea7644f004e0ab76cc99501dd47[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Oct 1 17:35:54 2014 -0400

    updated visualization scripts

[33mcommit 293275fabf6d52a33d2dbf93917b45f0cb0017c1[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Oct 1 13:01:50 2014 -0400

    fred_update script
    
    Fixed output format.

[33mcommit 56849e9ab578c09aa67d203d3a99b9ad2d9a40b4[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Oct 1 11:55:48 2014 -0400

    fred_csv updated
    
    Improved readability in spreadsheet data.

[33mcommit b85c6917878e08f32c75372a6501ccdbc5ed10df[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Oct 1 11:32:43 2014 -0400

    Remove excluded file

[33mcommit d5b66feecaef0995d0a3d4d79fb135258cbf4b55[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Oct 1 10:22:54 2014 -0400

    Makefile updated
    
    Compile single-threaded by default.

[33mcommit c9ea5a6a73f96272af2163c5e1ad9c52e45c69b3[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Oct 1 09:49:27 2014 -0400

    doc Makefile updated

[33mcommit 3703b42bd04e108e29fb326b96925988fb7354f9[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Wed Oct 1 09:45:07 2014 -0400

    Added fred_update script.
    
    Assuming that FRED was installed via git, this script pulls the current
    version of the original repository.

[33mcommit 68a91846d34ea22ce8da5cf8c520dc085a3b14c9[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Sep 30 17:11:43 2014 -0400

    Fixed round-off error bug.

[33mcommit a30334b46019db8820e48a185273abb5afb03c05[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Sep 30 17:06:42 2014 -0400

    RESULTS directory
    
    1. provide an initial RESULTS directory.
    2. improved handling of empty RESULTS directory.

[33mcommit c26817f0b6e1c81c4ac8d558f10d2f4bde7d54bb[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Sep 30 16:56:19 2014 -0400

    Removed obsolete file src/do_states

[33mcommit 154060c7aab27bbe040aef54f0cadf537cee0233[m
Author: JohnGrefenstette <gref@pitt.edu>
Date:   Tue Sep 30 16:55:11 2014 -0400

    Removed obsolete tutorial files.

[33mcommit 78678767da4da9b2a9c3c7a8d75e1d59dcb975f7[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Sep 30 09:51:03 2014 -0400

    compare.out files
    
    Created files that contained the output of a successful regression
    test, for comparison to later test outputs.

[33mcommit ea3d46e35af70e2b19ead0578b205bfbdcea9a14[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Sep 30 09:47:50 2014 -0400

    rt command updated
    
    rt command now gives clear message about tests passing or failing.

[33mcommit 2adf117596f994d8c62f0183b923890f1482d4a4[m
Author: John Grefenstette <gref@pitt.edu>
Date:   Tue Sep 30 09:24:16 2014 -0400

    Tutorial exercises
    
    Exercises to run FRED, and analyze output.

[33mcommit 0c4171299e316491540537aced31e6d7f256b0e3[m
Author: fred-publichealth-pitt-edu <fred.publichealth.pitt.edu@gmail.com>
Date:   Mon Sep 29 14:40:42 2014 -0400

    Update README.md

[33mcommit d32d77c4eba3079715e5167faccc911be8f4f6ac[m
Author: fred-publichealth-pitt-edu <fred.publichealth.pitt.edu@gmail.com>
Date:   Mon Sep 29 14:40:17 2014 -0400

    Update README.md

[33mcommit be5fa9073111502e1dc84f0d277b7fe7b7739da8[m
Author: fred-publichealth-pitt-edu <fred.publichealth.pitt.edu@gmail.com>
Date:   Mon Sep 29 14:37:59 2014 -0400

    Update README.md

[33mcommit 95118772eb104036a48adb28e3f26f13314ea470[m
Author: fred-publichealth-pitt-edu <fred.publichealth.pitt.edu@gmail.com>
Date:   Mon Sep 29 14:37:31 2014 -0400

    Update README.md

[33mcommit 399fd502e8780c3a05113cfca91a9e2445b916b7[m
Author: fred-publichealth-pitt-edu <fred.publichealth.pitt.edu@gmail.com>
Date:   Mon Sep 29 14:33:04 2014 -0400

    Update README.md

[33mcommit 89240e1a2b7ac8f3a650e2ad5590a1e2711d6913[m
Author: fred-publichealth-pitt-edu <fred.publichealth.pitt.edu@gmail.com>
Date:   Mon Sep 29 14:31:21 2014 -0400

    Update README.md

[33mcommit 361fa35eb8aa253a261979606310e65b3e8cd677[m
Author: fred-publichealth-pitt-edu <fred.publichealth.pitt.edu@gmail.com>
Date:   Mon Sep 29 14:14:59 2014 -0400

    Update README.md

[33mcommit 66b13413431a08504b2d534d0ace9a681bdb0efb[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Mon Sep 29 13:09:50 2014 -0400

    code optimization of FRED scripts

[33mcommit 4808285ac782018ee2813f72f7b98727bde121e7[m
Author: David Galloway <gallowdd@hotmail.com>
Date:   Mon Sep 29 08:55:17 2014 -0400

    Copy CVS (29-Sep-2014) to Github
    
    This is the first push of the code from CVS (29-Sep-2014) to the new
    Github repository

[33mcommit cab60e53bffba411232989472f93771213ad2259[m
Author: David Galloway <ddg5@phdl5.gsph.pitt.edu>
Date:   Mon Sep 29 07:38:33 2014 -0400

    first commit
