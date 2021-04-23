<?php
// Demo of ws_rbinomial and ws_binomial -- generating binomial distributed "successs" or "rates"

require('ws_binomial.php' );    // or whereer you chose to store the  wSurvey\Math release.
require('ws_rbinomial.php' );

use \wSurvey\math\binomial as lBin ;
use \wSurvey\math\rbinomial as lRbin ;

// initialize a binomial with N=1000 and rate=0.05
$N=1000;
$rate=0.05;
$ndo=40  ; // number of random variables to to draw
$binObj=lBin\binomInit($N,$rate);    // initialize a "binomial distribtion" array


$ameanRv=array_sum($aa)/$ndo   ;      // some stats for this binomial distribution
$amean=lBin\binomMean($binObj);
$asd=lBin\binomSd($binObj);
$askew=lBin\binom_skewness($binObj);
$akurt=lBin\binom_kurtosis($binObj);

$res1=[];                           // draw binomial random variables (counts of successes)
for ($ido=0;$ido<$ndo;$ido++) {
  $aa[]=lBin\binomRand($binObj);
}

$mess="Simple demo of binomial, and reverse binomial, random variable generation<p> ";
$mess.="$ndo binomial draws of # of successes (rate=$rate, N=$N).";
$mess.="<br>Mean: $amean,  sd=$asd, skew=$askew, Kurtosis=$akurt.  Average of draws=$ameanRv : ";
$t1=implode(', ',$aa);
$mess.='<br><tt>'.$t1.'</tt><hr>';

//Now the reverse binomial. Draw rates (rates given successs and attempts
$N=1000;
$m=5;
$ndo=40  ; // number to draw
$rBinObj=lRbin\rbinom_rv_init($m,$N);   // intialize a "reverse binomial object" -- use the default specifications
$bb=lRbin\rbinom_rv($ndo,$rBinObj);   // draw reverse binomial rvs

$mess.="$ndo reverse binomial draws of rates (m=$m N=$N) -- sorted.<br> Summary: ".$bb['message'] ;
$rgot=$bb['rvs'];             // note that the list of random variables is in the 'rv' index.
sort($rgot,SORT_NUMERIC);     // just for display purposes....
$t1=implode(', ',$rgot);
$amean=array_sum($rgot)/$ndo;
$mess.='<br>Average of draws: '.$amean.'  Median= '.$rgot[intval($ndo/2)] ;
$mess.='<br><tt>'.$t1.'</tt><hr>';


echo $mess;


?>
