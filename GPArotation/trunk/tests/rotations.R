#   Tests here only compare against previously computed values with this code
#   to ensure there was no accidental change. It would be better to have
#   comparisons with known correct values.


 Sys.getenv("R_LIBS")
 library()
 require("GPArotation")
 search()
 Sys.info()

require("stats")  
require("GPArotation")  

fuzz <- 1e-6 
all.ok <- TRUE  


  data(ability.cov)
  L <- loadings(factanal(factors = 2, covmat=ability.cov))
  

 if( 0.001 < max(abs(varimax(L, normalize=FALSE)$loadings -
          Varimax(L, normalize=FALSE)$loadings))) {
    cat("Calculated difference exceeds tolerance\n")
    cat("difference:\n")
    print(varimax(L, normalize=FALSE)$loadings -
          Varimax(L, normalize=FALSE)$loadings, digits=18)
    all.ok <- FALSE  
    } 

 if( 0.01 < max(abs(varimax(L, normalize=TRUE)$loadings -
          Varimax(L, normalize=TRUE, eps=1e-5)$loadings))) {
    cat("Calculated difference exceeds tolerance\n")
    cat("difference:\n")
    print(varimax(L, normalize=TRUE)$loadings -
          Varimax(L, normalize=TRUE, eps=1e-5)$loadings, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
           0.3863615904740822504,  0.4745127741495974161,
          -0.0110059418769087539,  0.6458720769633764514,
          -0.0262926272350604423,  0.8961141105684561348,
          -0.0180200526810754824,  0.4882928281695405048,
           0.9900944939102318543, -0.0370718282544326011,
           0.7905657274265397438,  0.0526109550054999417
      ), 2, 6))
 

  v <- oblimin(L)$loadings 
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
           0.3863615904740822504,  0.4745127741495974161,
          -0.0110059418769087539,  0.6458720769633764514,
          -0.0262926272350604423,  0.8961141105684561348,
          -0.0180200526810754824,  0.4882928281695405048,
           0.9900944939102318543, -0.0370718282544326011,
           0.7905657274265397438,  0.0526109550054999417
      ), 2, 6))

  v <- quartimin(L)$loadings
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
  	  0.551529228817982942, 0.4905002767031292898,
  	  0.217748645523411000, 0.6027046291262584399,
  	  0.291173432863349457, 0.8348885228488550636,
  	  0.154994397662456290, 0.4544843569140373241,
  	  0.969702339393929247, 0.0850652965070581996,
  	  0.803390575440818822, 0.1448091121037717866
      ), 2, 6))

  v <- targetT(L, Target=matrix(c(rep(1,3),rep(0,6),rep(1,3)), 6,2),
               eps=1e-5)$loadings 
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
  	  0.735795682866631218, 0.565351705145453853,
  	  0.433590223819374398, 0.664644550038417159,
  	  0.589924557708411568, 0.920006940799857786,
  	  0.317543426981046928, 0.500590650032113116,
  	  1.021758247914384077, 0.155121528590726393,
  	  0.872521244896209747, 0.208735706420634437
      ), 2, 6))

  v <- targetQ(L, Target=matrix(c(rep(1,3),rep(0,6),rep(1,3)), 6,2),
               eps=1e-5)$loadings  
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
          0.37067889993474656407, 0.638257130653133720,
          0.01855112570739854416, 0.640564749523800270,
          0.01576132191496706567, 0.884065831441111172,
          0.00524531003824213384, 0.480158078874985073,
          0.89458633399812259590, 0.383762977265515448,
          0.71793428958051475064, 0.388556883222951677
      ), 2, 6))

  v <- pstT(L, W = matrix(c(rep(.4,6),rep(.6,6)), 6,2),
           Target= matrix(c(rep(1,3),rep(0,6),rep(1,3)), 6,2),
               maxit=1000, eps=1e-5)$loadings     
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
          0.573125161748393785, 0.700868331877288475,
          0.214899397066479453, 0.681727425525818886,
          0.286558275327103040, 0.940272379393286339,
          0.152257795885557295, 0.510481967637567036,
          1.029289798076480578, 0.462598702071116141,
          0.850691132520651205, 0.456859727346562328
      ), 2, 6))

  v <- pstQ(L, W = matrix(c(rep(.4,6),rep(.6,6)), 6,2),
           Target= matrix(c(rep(1,3),rep(0,6),rep(1,3)), 6,2),
               maxit=1000, eps=1e-5)$loadings     
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
           -8111059.94622692652,  8111060.62253121007,
            1495036.43465861562, -1495035.79614594672,
            2331634.63904705830, -2331633.75893370388,
            1356735.91680212389, -1356735.43916810025,
          -23187491.19758165255, 23187491.68068471923,
          -18357040.58573083207, 18357041.05348757654
      ), 2, 6))

  v <- oblimax(L)$loadings  
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
          -0.595241134680514428, -0.593413784323445381,
          -0.648726078098998959, -0.257389315895841408,
          -0.896758845688784056, -0.345816215796869053,
          -0.487533402896166868, -0.184649783402633882,
          -0.260349007963903922, -0.999438062373664504,
          -0.290878679630171777, -0.832224919261500662
      ), 2, 6))

  v <- entropy(L, maxit=1000, eps=1e-5)$loadings  
 # this still does not converge
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
  	  0.534714740804540178, 0.508778102568043678,
  	  0.197348140750149392, 0.609689309353509956,
  	  0.262919828098457153, 0.844212045390758559,
  	  0.139616102327241837, 0.459441658926639795,
  	  0.966291466215733252, 0.117641548844535412,
  	  0.798063848020893585, 0.171756193883937508
      ), 2, 6))

  v <- quartimax(L, eps=1e-5)$loadings
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
          0.515866523962843160, 0.527879475961036904,
          0.175054634278874244, 0.616460231981747930,
          0.232057748479543163, 0.853211588623112749,
          0.122822468397975171, 0.464213243286899446,
          0.961376376417989453, 0.152689863976982837,
          0.791292800869773050, 0.200653429940987366
      ), 2, 6))

  v <- Varimax(L)$loadings  
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
  	   0.3384175759313114429, 0.508414890494446547464,
  	  -0.0654601124161610648, 0.670992229004664153535,
  	  -0.1016231721735353366, 0.930535379393095940515,
  	  -0.0589933707274080121, 0.506904360351960181497,
  	   0.9733094402675376289, 0.000234046050254643859,
  	   0.7702037184085044341, 0.085651123319384916965
      ), 2, 6))

  v <- simplimax(L, eps=1e-5)$loadings
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
          0.523583611303327312, 0.520226117818945788,
          0.184113022124463677, 0.613815719643687197,
          0.244596116053327067, 0.849702038129718673,
          0.129644684715025493, 0.462354355134084738,
          0.963520501269179652, 0.138517057902201340,
          0.794161628656258278, 0.188979901644201559
      ), 2, 6))

  v <- bentlerT(L)$loadings 
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
           0.3801726240258240241,  0.4741208368044214638,
          -0.0223632969057368826,  0.6514196922540864687,
          -0.0421105927111659756,  0.9039359851665277334,
          -0.0266594447192576613,  0.4925968005718689424,
           0.9961524457620027917, -0.0485973498906049697,
           0.7939648477384558811,  0.0440983921679098251
      ), 2, 6))

  v <- bentlerQ(L)$loadings 
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
          0.842839618436879490, 1.147474221357566826,
          0.966496634451572856, 0.899179738613143820,
          1.337254014063740337, 1.233960563362769269,
          0.727426977709420441, 0.667827007601493383,
          0.277572186047747749, 1.169666309336891752,
          0.345129603963082010, 1.050133434159491896
      ), 2, 6))

  v <- tandemI(L, maxit=1000, eps=1e-5)$loadings  
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
          -19.9573339993058916, -19.9414309667870207,
          -15.2547723459692186, -15.5482769465485902,
          -20.9207010736416947, -21.3343935951843306,
          -11.3177711730793167, -11.5453602088763034,
          -21.0732738925002536, -20.4756229416411912,
          -18.7993114151414424, -18.3586393947653654
      ), 2, 6))

  v <- tandemII(L, maxit=1000, eps=1e-5)$loadings 
 # this still does not converge
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
  	  0.572197044101002361, 0.4662247895688098054,
  	  0.243573415560656120, 0.5927388411683653935,
  	  0.326956608263186954, 0.8215352639437966120,
  	  0.174476792179181994, 0.4473668997335142894,
  	  0.972471249855535680, 0.0431091626026945812,
  	  0.808894688433769660, 0.1099794466209375043
      ), 2, 6))

  v <- geominT(L, eps=1e-5)$loadings  
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
           0.39672053553904490508,  0.4713295988080449250,
           0.00424452688463150020,  0.6389466007374070555,
          -0.00510976786312981532,  0.8864521406378518265,
          -0.00646959173137159373,  0.4830101828530461994,
           0.98709860078485589518, -0.0318959930081098297,
           0.79011178369962709045,  0.0558689642678330683
      ), 2, 6))

  v <- geominQ(L, eps=1e-5)$loadings  
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
          0.534721263659975854, 0.508771247100584523,
          0.197355957387199576, 0.609686779159006154,
          0.262930651479430233, 0.844208674501022327,
          0.139621992686633722, 0.459439868910532512,
          0.966292974385164483, 0.117629160286744874,
          0.798066049992627313, 0.171745962120156664
      ), 2, 6))

  v <- cfT(L)$loadings	
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
           0.3863615904740822504,  0.4745127741495974161,
          -0.0110059418769087539,  0.6458720769633764514,
          -0.0262926272350604423,  0.8961141105684561348,
          -0.0180200526810754824,  0.4882928281695405048,
           0.9900944939102318543, -0.0370718282544326011,
           0.7905657274265397438,  0.0526109550054999417
      ), 2, 6))

  v <- cfQ(L)$loadings	
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
  	  0.495330443338021176, 0.547195361446864537,
  	  0.151384273205308784, 0.622695868320644275,
  	  0.199304253086364791, 0.861451466010626055,
  	  0.105004533733904976, 0.468565194910632365,
  	  0.954843809781045660, 0.189293503899924942,
  	  0.783052579543945471, 0.230726576980168713
      ), 2, 6))

  v <- infomaxT(L, eps=1e-5)$loadings 
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
  	   0.39327554287862442894,  0.4693137508305071925,
  	  -0.00319802321222481794,  0.6422985517185823001,
  	  -0.01549245038490981718,  0.8912279460026399924,
  	  -0.01214605901641467763,  0.4856544522916727002,
  	   0.99260028929193111491, -0.0433225495465055510,
  	   0.79356458059567791530,  0.0471559021503157039
      ), 2, 6))

  v <- infomaxQ(L, eps=1e-5)$loadings 
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


  tst <- t(matrix(c(
          0.499563849536320603, 0.543711314635863419,
          0.156214100300119313, 0.621620013394660775,
          0.205986246705249010, 0.860034043856927166,
          0.108639140267620352, 0.467818061529905960,
          0.956284901268870557, 0.182598340120454394,
          0.784820467897764273, 0.225234145534117969
      ), 2, 6))

  v <- mccammon(L, maxit=1000, eps=1e-5)$loadings 
 # this still does not converge
  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test x. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 


 
cat("tests completed.\n")


if (! all.ok) stop("some tests FAILED")
