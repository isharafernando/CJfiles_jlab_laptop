CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION OFF_fmKP_F1 (x,wfn,dRN)
C   Polynomial fit to calculated off-shell correction in fmKP model
C   These off-shell corrections were calculated in a generalized
C   model with seperate valence, sea, and gluon corrections
C
C  Defined such that F1d = F1d(conv) + del^off F1d 
C       with OFF_fmKP = del^off F1d / F1d
C                    = off-shell correction w.r.t. F1d
C
C   Code and fits by L. Brady (August 2013)
C   Code automatically generated in python
C   email Lucas_Brady@hmc.edu if you want the python script to generate this
C  wfn: 
C       1 (Paris)
C       2 (AV18)
C       3 (CDBonn)
C       4 (WJC1)
C       5 (WJC2)
C  dRN:    [% change in nucleon radius in the deuteron]
C       1 (0.3%)
C       2 (0.6%)
C       3 (0.9%)
C       4 (1.2%)
C       5 (1.5%)
C       6 (1.8%)
C       7 (2.1%)
C       8 (2.4%)
C       9 (2.7%)
C       10 (3.0%)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Implicit None
      double precision OFF_fmKP_F1, x
      integer wfn, dRN
      double precision p(0:10)

      OFF_fmKP_F10.D0
      if (wfn.LT.1 .OR. wfn.GT.4) return
      if (dRN.LT.1 .OR. dRN.GT.10) return
!.......................................................................
      if (wfn.EQ.1) then     ! Paris
         if (dRN.EQ.1) then  ! 0.3%
            p(0) = 27.0128825948
            p(1) = -1613.26482632
            p(2) = 29369.7904368
            p(3) = -240904.431973
            p(4) = 1076004.74852
            p(5) = -2919745.3393
            p(6) = 5048922.25153
            p(7) = -5604183.92969
            p(8) = 3867520.51918
            p(9) = -1509842.57
            p(10) = 254336.474961
         else if (dRN.EQ.2) then  ! 0.6%
            p(0) = 25.0055108745
            p(1) = -1656.10699446
            p(2) = 31149.3543272
            p(3) = -259041.544177
            p(4) = 1156497.4648
            p(5) = -3110517.73009
            p(6) = 5294481.98587
            p(7) = -5743510.41889
            p(8) = 3841581.5827
            p(9) = -1437982.16236
            p(10) = 228761.452107
         else if (dRN.EQ.3) then  ! 0.9%
            p(0) = 18.609839817
            p(1) = -1493.8506941
            p(2) = 29435.3074529
            p(3) = -246937.255414
            p(4) = 1084348.44402
            p(5) = -2820489.34483
            p(6) = 4566244.50531
            p(7) = -4614406.87754
            p(8) = 2788626.15102
            p(9) = -895701.190083
            p(10) = 110028.683319
         else if (dRN.EQ.4) then  ! 1.2%
            p(0) = 6.7250528979
            p(1) = -1063.43273543
            p(2) = 23025.9710792
            p(3) = -193418.549291
            p(4) = 800310.881218
            p(5) = -1856273.11078
            p(6) = 2462115.07991
            p(7) = -1682876.98684
            p(8) = 270179.593497
            p(9) = 319795.176649
            p(10) = -142259.92713
         else if (dRN.EQ.5) then  ! 1.5%
            p(0) = -13.9045267797
            p(1) = -216.734403316
            p(2) = 9437.0570991
            p(3) = -77191.924148
            p(4) = 197633.98829
            p(5) = 116668.461373
            p(6) = -1692462.97145
            p(7) = 3925829.94622
            p(8) = -4418541.22357
            p(9) = 2529473.35518
            p(10) = -591225.989228
         else if (dRN.EQ.6) then  ! 1.8%
            p(0) = -52.00048867
            p(1) = 1443.65580327
            p(2) = -17939.6915624
            p(3) = 157574.777999
            p(4) = -998611.384184
            p(5) = 3942781.48966
            p(6) = -9564322.97386
            p(7) = 14326399.2514
            p(8) = -12944118.1874
            p(9) = 6475804.42528
            p(10) = -1379744.96852
         else if (dRN.EQ.7) then  ! 2.1%
            p(0) = -112.013321911
            p(1) = 4179.31251402
            p(2) = -64014.0630566
            p(3) = 555541.059974
            p(4) = -3020133.32078
            p(5) = 10363414.4301
            p(6) = -22667450.3255
            p(7) = 31497255.4785
            p(8) = -26907456.969
            p(9) = 12889305.674
            p(10) = -2651521.07315
         else if (dRN.EQ.8) then  ! 2.4%
            p(0) = 111.062906374
            p(1) = -4760.20435877
            p(2) = 65295.1727824
            p(3) = -387922.973833
            p(4) = 954978.616812
            p(5) = 20494.3084954
            p(6) = -5595698.81465
            p(7) = 13711816.96
            p(8) = -15699052.015
            p(9) = 9056405.69748
            p(10) = -2122902.58142
         else if (dRN.EQ.9) then  ! 2.7%
            p(0) = -378.408075561
            p(1) = 16381.246891
            p(2) = -269288.791948
            p(3) = 2318416.75533
            p(4) = -11882867.5498
            p(5) = 38194034.5083
            p(6) = -78835221.5896
            p(7) = 104330265.275
            p(8) = -85547565.1814
            p(9) = 39567181.1979
            p(10) = -7892491.1196
         else if (dRN.EQ.10) then  ! 3.0%
            p(0) = 54.2740967504
            p(1) = -3.37768929919
            p(2) = -37235.6919177
            p(3) = 619640.137597
            p(4) = -4557404.09759
            p(5) = 18351696.257
            p(6) = -44164271.3578
            p(7) = 65360963.1203
            p(8) = -58396529.6622
            p(9) = 28929461.5705
            p(10) = -6108217.58319
         endif
!.......................................................................
      else if (wfn.EQ.2) then     ! AV18
         if (dRN.EQ.1) then  ! 0.3%
            p(0) = 27.3591229679
            p(1) = -1635.43552006
            p(2) = 29734.6537385
            p(3) = -243714.490175
            p(4) = 1087799.39715
            p(5) = -2949813.27201
            p(6) = 5097424.59707
            p(7) = -5653495.39446
            p(8) = 3897507.25775
            p(9) = -1519348.7038
            p(10) = 255391.082087
         else if (dRN.EQ.2) then  ! 0.6%
            p(0) = 24.889495546
            p(1) = -1660.17425457
            p(2) = 31188.7290593
            p(3) = -258831.161103
            p(4) = 1151658.70093
            p(5) = -3083583.61012
            p(6) = 5218350.75093
            p(7) = -5618885.68777
            p(8) = 3721552.59069
            p(9) = -1374575.06679
            p(10) = 214521.779091
         else if (dRN.EQ.3) then  ! 0.9%
            p(0) = 17.0595341823
            p(1) = -1428.70908312
            p(2) = 28220.9412034
            p(3) = -235338.274131
            p(4) = 1019429.32666
            p(5) = -2596803.66043
            p(6) = 4077810.59064
            p(7) = -3936276.7127
            p(8) = 2208027.72898
            p(9) = -615913.748705
            p(10) = 51885.0140726
         else if (dRN.EQ.4) then  ! 1.2%
            p(0) = 1.52406120079
            p(1) = -822.219400912
            p(2) = 18685.1164979
            p(3) = -153876.739985
            p(4) = 590505.591374
            p(5) = -1165930.90518
            p(6) = 1010920.22727
            p(7) = 268831.377555
            p(8) = -1355249.54394
            p(9) = 1083662.64752
            p(10) = -297250.514353
         else if (dRN.EQ.5) then  ! 1.5%
            p(0) = 1.84299867418
            p(1) = -778.914898788
            p(2) = 16327.1968747
            p(3) = -117168.705997
            p(4) = 313153.551243
            p(5) = -11819.1410393
            p(6) = -1847094.22682
            p(7) = 4601067.33023
            p(8) = -5304751.4026
            p(9) = 3074386.22863
            p(10) = -724026.074628
         else if (dRN.EQ.6) then  ! 1.8%
            p(0) = -88.2515278376
            p(1) = 3107.57864778
            p(2) = -46222.8106963
            p(3) = 402168.301989
            p(4) = -2233801.82979
            p(5) = 7834915.49319
            p(6) = -17444952.4785
            p(7) = 24580592.3355
            p(8) = -21232425.349
            p(9) = 10263667.3263
            p(10) = -2127875.20561
         else if (dRN.EQ.7) then  ! 2.1%
            p(0) = -179.73041304
            p(1) = 7262.39142628
            p(2) = -116082.336335
            p(3) = 1004306.71518
            p(4) = -5281219.76373
            p(5) = 17477408.4259
            p(6) = -37055912.7384
            p(7) = 50201846.7098
            p(8) = -42010987.9824
            p(9) = 19783346.5374
            p(10) = -4010958.72489
         else if (dRN.EQ.8) then  ! 2.4%
            p(0) = -350.475189493
            p(1) = 15172.7428178
            p(2) = -249641.955652
            p(3) = 2150977.39395
            p(4) = -11031605.8342
            p(5) = 35474242.5942
            p(6) = -73246854.1931
            p(7) = 96963190.1531
            p(8) = -79529875.0149
            p(9) = 36795803.1773
            p(10) = -7342542.69434
         else if (dRN.EQ.9) then  ! 2.7%
            p(0) = -359.879659421
            p(1) = 15577.7076267
            p(2) = -257988.735334
            p(3) = 2244817.34238
            p(4) = -11645439.5257
            p(5) = 37885003.3489
            p(6) = -79112919.4471
            p(7) = 105876551.641
            p(8) = -87760286.7007
            p(9) = 41020503.5107
            p(10) = -8267291.0297
         else if (dRN.EQ.10) then  ! 3.0%
            p(0) = -321.813599288
            p(1) = 13996.6668787
            p(2) = -234705.200031
            p(3) = 2076769.76813
            p(4) = -10979486.1926
            p(5) = 36368590.3124
            p(6) = -77210642.9469
            p(7) = 104905546.204
            p(8) = -88184010.0719
            p(9) = 41766254.4298
            p(10) = -8524204.58489
         endif
!.......................................................................
      else if (wfn.EQ.3) then     ! CDBonn
         if (dRN.EQ.1) then  ! 0.3%
            p(0) = 24.8520401588
            p(1) = -1525.63486852
            p(2) = 28181.4080417
            p(3) = -232655.503982
            p(4) = 1042270.6971
            p(5) = -2831561.32794
            p(6) = 4896629.73336
            p(7) = -5430994.8251
            p(8) = 3743157.20151
            p(9) = -1459064.43987
            p(10) = 245457.502472
         else if (dRN.EQ.2) then  ! 0.6%
            p(0) = 23.0310692895
            p(1) = -1591.13260067
            p(2) = 30583.3290434
            p(3) = -257018.802337
            p(4) = 1155030.98912
            p(5) = -3122916.59196
            p(6) = 5342117.03501
            p(7) = -5826595.05336
            p(8) = 3923022.12298
            p(9) = -1481691.40187
            p(10) = 238875.402401
         else if (dRN.EQ.3) then  ! 0.9%
            p(0) = 18.5112545023
            p(1) = -1530.73621968
            p(2) = 30844.5931405
            p(3) = -262888.209107
            p(4) = 1174657.43912
            p(5) = -3121685.00438
            p(6) = 5196686.87417
            p(7) = -5454742.33107
            p(8) = 3483336.45322
            p(9) = -1221370.7437
            p(10) = 176423.015344
         else if (dRN.EQ.4) then  ! 1.2%
            p(0) = 10.4708476416
            p(1) = -1311.53197698
            p(2) = 28470.7311342
            p(3) = -246331.60094
            p(4) = 1082235.29823
            p(5) = -2769468.1016
            p(6) = 4342143.50848
            p(7) = -4159703.16079
            p(8) = 2295664.29232
            p(9) = -617981.271163
            p(10) = 45922.498417
         else if (dRN.EQ.5) then  ! 1.5%
            p(0) = -1.66135308362
            p(1) = -884.703822038
            p(2) = 22367.2596206
            p(3) = -196386.573096
            p(4) = 817707.051427
            p(5) = -1868198.43728
            p(6) = 2367651.45802
            p(7) = -1401283.63766
            p(8) = -77157.2395277
            p(9) = 527197.391498
            p(10) = -191469.831946
         else if (dRN.EQ.6) then  ! 1.8%
            p(0) = -21.5373424344
            p(1) = -80.3697142981
            p(2) = 9778.07038366
            p(3) = -90524.9772073
            p(4) = 273576.400825
            p(5) = -96551.51166
            p(6) = -1346923.04668
            p(7) = 3593048.37251
            p(8) = -4235175.3019
            p(9) = 2478463.56846
            p(10) = -586172.014386
         else if (dRN.EQ.7) then  ! 2.1%
            p(0) = -37.6298323842
            p(1) = 502.773692481
            p(2) = 851.309368225
            p(3) = -12407.6347772
            p(4) = -156536.370649
            p(5) = 1403484.92272
            p(6) = -4687229.52836
            p(7) = 8319105.87214
            p(8) = -8342607.954
            p(9) = 4477382.93268
            p(10) = -1003230.20522
         else if (dRN.EQ.8) then  ! 2.4%
            p(0) = -16.8197288599
            p(1) = 116.265829402
            p(2) = 490.551591898
            p(3) = 36617.6677798
            p(4) = -595067.426675
            p(5) = 3274581.1387
            p(6) = -9282402.30549
            p(7) = 15160384.8828
            p(8) = -14450803.2451
            p(9) = 7490704.25354
            p(10) = -1635484.89588
         else if (dRN.EQ.9) then  ! 2.7%
            p(0) = -160.181138499
            p(1) = 5939.23058354
            p(2) = -88498.8074302
            p(3) = 745189.783919
            p(4) = -3951576.5457
            p(5) = 13323794.6377
            p(6) = -28788802.7836
            p(7) = 39649298.1481
            p(8) = -33638306.1531
            p(9) = 16020445.2617
            p(10) = -3278384.57804
         else if (dRN.EQ.10) then  ! 3.0%
            p(0) = -247.719310261
            p(1) = 9954.24208409
            p(2) = -156236.8311
            p(3) = 1330158.09169
            p(4) = -6914506.55571
            p(5) = 22696802.6098
            p(6) = -47831535.5996
            p(7) = 64486281.6717
            p(8) = -53737095.4374
            p(9) = 25204999.3079
            p(10) = -5089846.43245
         endif
!.......................................................................
      else if (wfn.EQ.4) then     ! WJC1
         if (dRN.EQ.1) then  ! 0.3%
            p(0) = 35.3685582504
            p(1) = -2013.34445437
            p(2) = 35835.9303794
            p(3) = -292323.194838
            p(4) = 1311469.24413
            p(5) = -3598713.20151
            p(6) = 6327889.80649
            p(7) = -7178638.60389
            p(8) = 5088804.26306
            p(9) = -2051389.76266
            p(10) = 358893.981366
         else if (dRN.EQ.2) then  ! 0.6%
            p(0) = 37.9028291532
            p(1) = -2254.82219794
            p(2) = 40689.1872691
            p(3) = -335990.26801
            p(4) = 1520040.79687
            p(5) = -4195992.29078
            p(6) = 7406118.6947
            p(7) = -8414500.45511
            p(8) = 5958562.07164
            p(9) = -2392169.71644
            p(10) = 415186.537053
         else if (dRN.EQ.3) then  ! 0.9%
            p(0) = 37.0854254367
            p(1) = -2334.99061936
            p(2) = 42742.5906478
            p(3) = -355060.411247
            p(4) = 1602816.69747
            p(5) = -4392028.90002
            p(6) = 7661505.66999
            p(7) = -8564928.69131
            p(8) = 5937945.55673
            p(9) = -2319539.88687
            p(10) = 388431.089
         else if (dRN.EQ.4) then  ! 1.2%
            p(0) = 30.1020772397
            p(1) = -2124.38629769
            p(2) = 39788.6516086
            p(3) = -330303.131267
            p(4) = 1461992.52006
            p(5) = -3876505.25754
            p(6) = 6461620.41594
            p(7) = -6802098.41862
            p(8) = 4354352.82526
            p(9) = -1524799.49885
            p(10) = 217466.467611
         else if (dRN.EQ.5) then  ! 1.5%
            p(0) = 10.1584363908
            p(1) = -1304.80360851
            p(2) = 26466.529893
            p(3) = -216044.908627
            p(4) = 871066.556472
            p(5) = -1948796.33636
            p(6) = 2413097.69147
            p(7) = -1343890.22396
            p(8) = -208997.031543
            p(9) = 629133.769211
            p(10) = -221516.733433
         else if (dRN.EQ.6) then  ! 1.8%
            p(0) = -24.9730189175
            p(1) = 269.918230339
            p(2) = -279.086955765
            p(3) = 17984.0190179
            p(4) = -337918.408801
            p(5) = 1958158.68706
            p(6) = -5694291.51017
            p(7) = 9450737.1354
            p(8) = -9122751.7106
            p(9) = 4785272.39675
            p(10) = -1058166.35097
         else if (dRN.EQ.7) then  ! 2.1%
            p(0) = -102.75691584
            p(1) = 3810.77359842
            p(2) = -60016.8925782
            p(3) = 533616.736395
            p(4) = -2950052.77612
            p(5) = 10231973.617
            p(6) = -22544065.9277
            p(7) = 31502154.2272
            p(8) = -27044336.6347
            p(9) = 13017381.6335
            p(10) = -2691658.24872
         else if (dRN.EQ.8) then  ! 2.4%
            p(0) = -154.671631162
            p(1) = 6219.43621962
            p(2) = -101908.456315
            p(3) = 907197.161927
            p(4) = -4906562.00049
            p(5) = 16629656.0261
            p(6) = -35970271.4861
            p(7) = 49578340.1729
            p(8) = -42135984.9658
            p(9) = 20130264.0252
            p(10) = -4138432.71124
         else if (dRN.EQ.9) then  ! 2.7%
            p(0) = -186.25753541
            p(1) = 7641.75991112
            p(2) = -126758.094276
            p(3) = 1133758.18294
            p(4) = -6133243.46316
            p(5) = 20790237.8147
            p(6) = -45028236.3265
            p(7) = 62217281.3018
            p(8) = -53058195.0979
            p(9) = 25451642.4568
            p(10) = -5255975.99393
         else if (dRN.EQ.10) then  ! 3.0%
            p(0) = -188.128167275
            p(1) = 7740.88243616
            p(2) = -129928.013705
            p(3) = 1178775.13469
            p(4) = -6471143.76818
            p(5) = 22234165.1291
            p(6) = -48753427.862
            p(7) = 68140768.0465
            p(8) = -58742669.928
            p(9) = 28472732.0308
            p(10) = -5939304.57853
         endif
!.......................................................................
      else if (wfn.EQ.5) then     ! WJC2
         if (dRN.EQ.1) then  ! 0.3%
            p(0) = 29.6286245116
            p(1) = -1738.25866758
            p(2) = 31380.7898811
            p(3) = -256835.561613
            p(4) = 1148653.74002
            p(5) = -3127984.29396
            p(6) = 5438022.38004
            p(7) = -6078332.31682
            p(8) = 4230854.35186
            p(9) = -1668675.81391
            p(10) = 284499.347337
         else if (dRN.EQ.2) then  ! 0.6%
            p(0) = 28.7383442574
            p(1) = -1827.96935738
            p(2) = 33836.7825858
            p(3) = -280344.473599
            p(4) = 1255483.55073
            p(5) = -3401232.8589
            p(6) = 5850722.48908
            p(7) = -6435335.1422
            p(8) = 4380405.98602
            p(9) = -1676359.11772
            p(10) = 274382.651143
         else if (dRN.EQ.3) then  ! 0.9%
            p(0) = 23.2674790088
            p(1) = -1700.21655487
            p(2) = 32554.5160017
            p(3) = -271279.747017
            p(4) = 1196999.1459
            p(5) = -3151215.24777
            p(6) = 5198786.40168
            p(7) = -5399547.72823
            p(8) = 3397478.63233
            p(9) = -1163152.70714
            p(10) = 160684.692726
         else if (dRN.EQ.4) then  ! 1.2%
            p(0) = 10.4783286259
            p(1) = -1222.00707797
            p(2) = 25210.9892891
            p(3) = -209217.413198
            p(4) = 869175.096507
            p(5) = -2048784.6537
            p(6) = 2814831.91874
            p(7) = -2103017.76538
            p(8) = 581951.876438
            p(9) = 189809.931085
            p(10) = -119267.066036
         else if (dRN.EQ.5) then  ! 1.5%
            p(0) = -10.7457077782
            p(1) = -301.266770618
            p(2) = 9807.75261461
            p(3) = -74300.999486
            p(4) = 162467.397136
            p(5) = 272957.083296
            p(6) = -2077636.29118
            p(7) = 4498529.87708
            p(8) = -4931501.01847
            p(9) = 2785411.4126
            p(10) = -646116.371704
         else if (dRN.EQ.6) then  ! 1.8%
            p(0) = -58.0862178277
            p(1) = 1759.86749679
            p(2) = -24191.0928434
            p(3) = 216732.15741
            p(4) = -1314298.90032
            p(5) = 4976088.91107
            p(6) = -11718701.1264
            p(7) = 17199929.956
            p(8) = -15320098.198
            p(9) = 7586102.7642
            p(10) = -1604165.62685
         else if (dRN.EQ.7) then  ! 2.1%
            p(0) = -139.199395178
            p(1) = 5468.57107566
            p(2) = -86682.5639209
            p(3) = 755123.522255
            p(4) = -4035479.55522
            p(5) = 13572393.3508
            p(6) = -29173520.7098
            p(7) = 39970316.922
            p(8) = -33763139.4266
            p(9) = 16027216.4116
            p(10) = -3272703.89299
         else if (dRN.EQ.8) then  ! 2.4%
            p(0) = -277.451125962
            p(1) = 11786.0663174
            p(2) = -193013.836274
            p(3) = 1669122.41125
            p(4) = -8635054.51983
            p(5) = 28031614.3836
            p(6) = -58389998.859
            p(7) = 77904289.1725
            p(8) = -64346856.8062
            p(9) = 29961229.397
            p(10) = -6014287.97271
         else if (dRN.EQ.9) then  ! 2.7%
            p(0) = -266.738363899
            p(1) = 11530.3282594
            p(2) = -193306.60705
            p(3) = 1712194.61966
            p(4) = -9059011.4624
            p(5) = 29991574.1702
            p(6) = -63553907.3165
            p(7) = 86095342.5307
            p(8) = -72098645.9182
            p(9) = 33998237.1214
            p(10) = -6905527.46243
         else if (dRN.EQ.10) then  ! 3.0%
            p(0) = -320.855472137
            p(1) = 13750.0744838
            p(2) = -227857.619994
            p(3) = 1995489.41917
            p(4) = -10461721.571
            p(5) = 34426967.8853
            p(6) = -72712420.1823
            p(7) = 98388507.0121
            p(8) = -82433263.9
            p(9) = 38939471.8515
            p(10) = -7930766.47026
         endif
      endif

      OFF_fmKP_F1=(p(0)*x+p(1)*x**2+p(2)*x**3+p(3)*x**4
     &                + p(4)*x**5+p(5)*x**6+p(6)*x**7
     &                + p(7)*x**8+p(8)*x**9+p(9)*x**10
     &                + p(10)*x**11+p(11)*x**12)
     &                * (1D0/(1D0+exp(-0.15D0*(x-50D0))))
      return
      end


