<?php
namespace wSurvey\Math\binomial ;

use wSurvey\Math\logGamma as lnG ;
use  wSurvey\Math\sortSearch  as ss;

// 22 APril 2021
// binomial distribution functions. This uses functions from the gburtini library
// $bObj=binomInit($n,$p)  : initialzie
//    n=# trials, p=prob success.
//    bObj is used in other functions
// aProb=binomPdf($bObj,$k,$puse)      : lookup pdf probablity
///   pdf of binomial. 2 or 3 argument version:
//   2 argument : bObj (as created by binomInit)-- the #trials and prob successes.
//              k : number of successes
//   3 arguments:
//      bObj : # of trials (integer)
//      k  : # of succese  (integer
//      puse : prob of success (float)
//
//   $rv=function binomRand($bObj,$useImp,$attempts,$nsd)     .... useImp, $attempts, and $nsd are optional
//     draw a single random binomial value
//      $bObj : as created by binomInit
//      $useImp: how to create
//        0 : use repeated bernoulli (slow, especifially if #trials is large)
//        1  use importance sampling. (faster)
//        2 : use importance sampling 2 steps -- first step uses CDF lookup in "fat" region (defined mode  +/=  m*sd)
//           Otherwise use importance sampling ON THE NON-FAT region. THIS IS TEH DEFAULT
//       3: choose 0 or 2 (if n<1000, use 0)
//    attempts: number of times to try if 1 or 2. If no success, use repeated bernoulli.  If 0, use max of 5 and  4*sqrt(n)
//    nsd:  number of sds (used in method 2) . < 0 means auto choose (typically about 3.1)
//


/* test ...
 $nTrials=2000000    ;
 $p=0.4 ;
 $ndo=69       ;
 $xx=binomInit($nTrials,$p);
 $amean=binomMean($xx);
 $avar=binomVariance($xx);
 $asd=binomSd($xx);
 echo "<br>mean=$amean  ,  sd=$asd  (variance=$avar) ";
 echo "<p> ";

// draw $ndo rvs
  $tim2a=microtime(true);
  $gooR=binomRandDrawImpN($ndo,$xx['p'],$xx['n']);
  $tim2b=microtime(true);
  $amean=array_sum($gooR)/count($gooR);
  $avar=compStdDev($gooR);
  $dt=$tim2b-$tim2a ;
  echo "<br>range Eps $dt mean=$amean var= $avar ";

// draw rvs one at a time
  $tim3a=microtime(true);
  for ($mm=0;$mm<$ndo;$mm++) {
     $r1=binomRand($xx,2,100,0);
     $goo[]=$r1;
  }
  $tim3b=microtime(true);
  $dt=$tim3b-$tim3a ;
  $amean=array_sum($goo)/count($goo);
  $avar=compStdDev($goo);
  echo "<br>range Eps $dt mean=$amean var=$avar ";
*/


//===========
// binomial distribution. non class based version (extenstion of gburtini package)

// initlize "object (associative array) used in other binomial calculations
// returns: [n, p, lng, errors]
//    errors is empty if no errors, otherwise an array of messages
//    lng=logGamma(n+1) (=log(n!))  -- useful to have this handy
// binomial distribution pdf= (n m) p^k  (1-p)^(n-k).
//  n: number of "events"
//  m: number of "successes" that occur during these n events
//  p : probablity of a single event yielding a success
// (n m) : choose m from n =  (n!) /(m!)(n-m)!
function  binomInit($n,$p) {
 $errs=[] ;
 $rets=['n'=>$n,'p'=>$p,'lng'=>null,'errs'=>$errs];
 if (is_numeric($n)) {
     $n=intVal($n);
     if ($n<0) {
         $errs[]='n (events) is lt 0 ';
     } else {
         $rets['lng']=lng\logGamma($n+1);
     }
  } else {
     $errs[]='n is non-numeric' ;
  }
 if (is_numeric($p)  ) {
     $p=floatVal($p);
     if ($p<0 || $p>1.0 ) $errs[]='p (success rates) is not between 0.0. and 1.0 ';
  } else {
     $errs[]='p is non-numeric' ;
  }
  $rets['errs']=$errs;

  return $rets;
}

// ----- pdf
function  binomPdf($xx,$k,$puse=null) {

  if (is_null($puse)) {
     $n = $xx['n'];
     $p = $xx['p'];
     $lng=$xx['lng'];
  } else {
      $n=$xx;
      $p=$puse;
      $lng=lng\logGamma($n+1);
   }
   $k=intval($k);
   $logBinom =  $lng- lng\logGamma($k+1)- lng\logGamma(1+$n-$k);
   $logP = $logBinom + $k * log($p) + ($n - $k)*log(1-$p);
   return exp($logP);
}

// ----- cdf (NOT efficient)
 function  binomCdf($xx,$k)   {
    $accumulator = 0.0;

    for ($i=0; $i<=$k; $i++) {
            $accumulator += binomPdf($xx,$i);
     }
    return $accumulator;
}

// inverse cdf (k such binomCdf(xx,k)=p  (NOT efficient)
function  binomIcdf($xx,$p)    {
  if ($p < 0 || $p > 1) {
       throw new InvalidArgumentException("Parameter (\$p = " . var_export($p, true) . " must be between 0 and 1. ");
   }
   if ($p == 1)    return INF;
   $accumulator = 0.0;
   $k = 0;
   do {
      $delta =binomPdf($xx,$k);
      $accumulator = $accumulator + $delta;
      if ($accumulator >= $p) return $k;
      $k = $k + 1;
   } while ($k < $this->n);
    return $k;
}

// -=------ mean, variance, sd

function  binomMean($xx)    {
     $n = $xx['n'];
     $p = $xx['p'];
    return $n * $p;
}
function  binomVariance($xx)    {
     $n = $xx['n'];
     $p = $xx['p'];
     return $n * $p * (1-$p);
}
function  binomSd($xx)    {
     $n = $xx['n'];
     $p = $xx['p'];
     $t1=binomVariance($xx);
     return sqrt($t1);
}

// -=------ skewness and kurtosis
function  binom_skewness($xx)  {
     $n = $xx['n'];
     $p = $xx['p'];
     $glug=sqrt($n * (1-$p) * $p) ;
     if ($glug==0) return 0 ;
     return (1-2 * $p)/$glug;
//     return (1-2 * $p)/sqrt($n * (1-$p) * $p);
}

function  binom_kurtosis($xx)   {
     $n = $xx['n'];
     $p = $xx['p'];
     if ($p==1.0 || $p==0.0) return 0;
     return 3+(1-6 * (1-$p) * $p)/($n * (1-$p) * $p);
}


//=========
// stupid (repeated bernouli) random
function binomRandSimple($n , $p) {
   $x = 0;
   for ($i = 0; $i < $n; $i++) {
     $gog=mt_rand()/mt_getrandmax();
     if ( $gog  < $p)   $x = $x + 1;
   }
   return $x  ;
}

//======
// draw a single random binomial value
// useIMP=0 : use repeated bernoulli
//        1  use importance sampling.
//        2 : use importance sampling 2 steps -- first step uses CDF lookup in "fat" region (defined mode  +/=  m*sd)
//           Otherwise use bernoulli ON THE NON-FAT region.
//       3: choose 0 or 2 (if n<1000, use 0)
// attempts: number of times to try if 1 or 2. If no success, use repeated bernoulli.  If 0, use max of 5 and  4*sqrt(n)
// nsd:  number of sds (used in method 2) . < 0 means auto choose

function binomRand($xx,$useImp=2,$attempts=0,$nsd=-1)    {
   $nUse = $xx['n'];
   $pUse = $xx['p'];

   if ($useImp==3) {        // auto deetermin
       if ($nUse <100 ) {      // dont bother with   importance if n is small
           $useImp=0;
       } else {
           $amean=$nUse * $pUse ;
           if ($nUse < 1000 && $amean>50) {   // if largish mean and not too big n, use bernoulli
              $useImp=0;
           } else {
              $useImp=2;
           }
        }
   }     // usimp=3

    if ($useImp==2 && $nsd<0) {   // auto determin # of sds
       if ($pUse>0.1) {
           if ($pUse>0.3) {
              $nsd=3.1;
           } else {
              $nsd=3.0;
           }
       } else {
          $nsd=3.05;
       }

    }      // determin nsd


 // draw random, using one of three tricks
      if ($useImp==0)   return binomRandSimple($nUse, $pUse);     // repeated bernoulli

      if ($useImp==1) {      // importance sampling
          $yval= binomRandDrawImp($nUse,  $pUse ,$attempts);
          if ($yval>=0) return $yval;
          return  binomRandSimple($nUse, $pUse);  // importance failed, use repeated bernoulli draws
      }

      if ($useImp==2) {  // importance sampling with fat area first
         $yval= binomRandDrawImp2($nUse, $pUse,0,$nsd) ; // ,$attempts,$nsd);
          if ($yval>=0) return $yval;              // if <0, abs(yval) = # of attempts
         $yval= binomRandDrawImp2($nUse, $pUse,0,$nsd*2) ; //   // expand the core
          if ($yval>=0) return $yval;              // if <0, abs(yval) = # of attempts
          return binomRandSimple($nUse,  $pUse );  // importance failed, use repeated bernoulli draws
      }

}


// importance sampling, using up to $attempts draws. If nothing found, return -$attempts
function binomRandDrawImp($n, $p,$attempts=0 )    {
    if ($p<=0) return 0;    // deal with simple cases
    if ($p>=1) return $n ;
    $maxRand=mt_getrandmax();
    if ($attempts==0)   $attempts=intVal(max(4*sqrt($n),5));  //  a guess as to where it might be better to give up and do the brutish force method

    $gammaN=lnG\logGamma($n + 1) ;
    $lnP=log($p); $lnP1=log(1-$p);

// the max pdf for this binomial.
     $yMode=floor(($n+1)*$p);

     $logBinom = $gammaN - lnG\logGamma($yMode + 1) -  lnG\logGamma(1+$n - $yMode);   // binomial pdf at mode (n*p)
     $logP = $logBinom + $yMode * $lnP + ($n - $yMode)*$lnP1;
     $maxPdf=exp($logP);

     for ($i = 0; $i < $attempts; $i++) {
          $x=mt_rand(0,$n);               // candidate value between 0 and n

          $logBinom = $gammaN - lnG\logGamma($x + 1) - lnG\logGamma(1+$n - $x );
           $logP = $logBinom + $x * $lnP  + ($n - $x)*$lnP1;
           $xpdf=exp($logP);

           $xpdfNorm=$xpdf/$maxPdf  ;   // normalise. note that maxPdf >= pdf of any possible $x
           $xtest=mt_rand()/$maxRand;
           if ($xtest<=$xpdfNorm)  return $x;
       }

       return -$attempts ;         // signals error (so use draw instead of drawImp
 }

// importance sampling, using split sample
 function binomRandDrawImp2($n, $p,$attempts=0,$nsds=1)    {
      if ($p<=0) return 0;    // deal with simple cases
      if ($p>=1) return $n ;

      if (!is_numeric($nsds) || $nsds<0) $nsds=1;  // ignore goofy stuff

      $maxRand=mt_getrandmax();
      if ($attempts==0)   $attempts=intVal(max(4*sqrt($n),5));          //  a guess as to where it might be better to give up and do the brutish force method

      $gammaN= lnG\logGamma ($n + 1) ;
      $lnP=log($p); $lnP1=log(1-$p);

// the max pdf for this binomial.
      $yMode=floor(($n+1)*$p);

// get cdf within nsds sds
      $sdx=sqrt($p*(1-$p)*$n)*$nsds;
      $Y0=max(0,round($yMode-$sdx));
      $Y1=min($n,round($yMode+$sdx));
      $icum=0;
      $cums=[];

    $logBinom = $gammaN - lnG\logGamma ($Y0 + 1) -  lnG\logGamma ($n - $Y0 + 1);
    $logP = $logBinom + $Y0 * $lnP  + ($n - $Y0)*$lnP1;
    $xpdf=exp($logP);
    $cand0=$xpdf ;
    $icum=$xpdf;
    $cums=[$icum];
    $pcomp=1-$p;

    for ($x=$Y0+1;$x<=$Y1;$x++) {
       $denom=($pcomp)*($x);
       $num=$p* ($n - ($x-1));
       $ff12=$num/$denom;
       $xpdf=$xpdf*$ff12;
        $icum+=$xpdf;
       $cums[]=$icum;
     }          // icum will have cdf of the entire dense region
     $dy=$Y1 - $Y0 ;

    $cand1=$xpdf;

//binarySearchIndexL: -1 : < cums[0], cums.length-1: > end(cums), jj : between cums[jj] and cums[jj+1]
// if -1: assign to cum[0], >0 then assign cums[jj+1], cums.length-1, jj+1 (should never happen)
// in dense region? do cdf match

    $asplit=mt_rand()/$maxRand;
    if ($asplit<=$icum)     {  // in dense region == use cdf lookup
       $afind=$icum*mt_rand()/$maxRand;
       $iat=ss\binarySearchMatrix($cums,$afind);
       if ($iat>=count($cums))   return $Y1; // should never happen
       if ($iat<0) return $Y0;
       return ($Y0+$iat+1);
    }

    $maxPdf=max($cand0,$cand1) ;  // max pdf of "non-fat region"

// not dense region. Use importance sampling (discard if falls in dense region)
    $iattempts=0  ;
    $t1=microtime(true);
    for ($i = 0; $i <1111; $i++) {
       if ($iattempts>$attempts) {
           return -$attempts ;  // give up
       }

       $x=mt_rand(0,$n);               // candidate value between 0 and n
       if ($x>=$Y0 && $x<=$Y1) continue ;  // avoicd double counting of "fat range" . This is kind of a stupid, but if fat region is relatively small, might be cheaper than a smarter method
       $iattempts++ ;

       $logBinom = $gammaN -  lnG\logGamma ($x + 1) -  lnG\logGamma ($n - $x + 1);
       $logP = $logBinom + $x * $lnP  + ($n - $x)*$lnP1;
       $xpdf=exp($logP);

       $xpdfNorm=$xpdf/$maxPdf  ;   // normalise. note that maxPdf >= pdf of any possible $x
       $xtest=mt_rand()/$maxRand;
       if ($xtest<=$xpdfNorm) {
         $t2=microtime(true);
         $dt=$t2-$t1 ;
           return $x;
       }
    }   // for

// should NEVER get here ($i> 111111111
    return -$attempts ;         // signals error (so use draw instead of drawImp
}


//==========
// imporatnce sampling using non uniform prior range. Return $ndo values 
// this is efficient if pulling multiple values from same binomial distribution (same ntrials and p)
// otherwise, it is slow (it does a fair amount of work at ffirst, but this work is used to quickly find rvs)
// #of rvs to return (in an array, probaility of success, # trials, resolution (of "trianglar cdf")
// $res sets the "range of x" -- small res means larger range (if binomPdf(x) < res, don't includein range)
//  setting this larger can speed things up a small amount, but will truncate range of values returned
function binomRandDrawImpN($ndo,$p,$nTrials,$res=1e-7) {
  $jj=1;
  $keeps=[];

  $asd=intVal(sqrt($p*(1-$p)*$nTrials));
  $xtry=intval($p*$nTrials);
  $xtry0a=$xtry-1;
  $xtry2a=$xtry+1;

  $p0=binomPdf($nTrials,$xtry0a,$p);
  $p1=binomPdf($nTrials,$xtry,$p);
  $p2=binomPdf($nTrials,$xtry2a,$p);

  $keeps[]=[$xtry0a,$p0,$xtry0a,$p0,$p0,$p0];
  $keeps[]=[$xtry,$p1,$xtry,$p1,$p1,$p1];
  $keeps[]=[$xtry2a,$p2,$xtry2a,$p2,$p2,$p2];

  $xtry=intVal($xtry0a-$asd);
  for (;;) {
    $jj++;
    if ($xtry<0) {
      $keeps[]=[0,0,0,0,0,0];
      break;
    }
    $xtry2=$xtry+$asd-1 ;
    $p1=binomPdf($nTrials,$xtry,$p);
    $p2=binomPdf($nTrials,$xtry2,$p);
    $mx=max($p1,$p2);
    $cc=$mx*$asd;
    $keeps[]=[$xtry,$p1,$xtry2,$p2,$mx,$cc];
    if ($p1<$res) break;
    $xtry=intVal($xtry-$asd) ;
  }    // low search

  $xtry=intVal($xtry2a+$asd);
  for (;;) {
    $jj++;
    if ($xtry>$nTrials) {
       $keeps[]=[$nTrials,0,0,0,0,0,0];
    }
    $xtry0=$xtry+1-$asd;
    $p0=binomPdf($nTrials,$xtry0,$p);
    $p1=binomPdf($nTrials,$xtry,$p);
    $mx=max($p0,$p1);
    $cc=$mx*$asd;
    $keeps[]=[$xtry0,$p0,$xtry,$p1,$mx,$cc];
    if ($p1<$res)  break;

    $xtry=intVal($xtry+$asd) ;
  }    // high search

// keeps : [0]:xLow value, [1] binom(xlow), [2] xhHigh, {3] binom(xhigh), 
//  [4] "raw bin pdf" :max of 1 and 3, [5] "raw bin cdf"= dx * max(1,3),
// [6] normalized bin pdf, [7] normalized bin cdf

  $xs=array_column($keeps,0);
  array_multisort($xs, SORT_ASC,SORT_NUMERIC,$keeps);   // lowest x to highest x

  $tot1=array_sum(array_column($keeps,5)) ;
  $cum1=0;
  for  ($ik=0;$ik<count($keeps);$ik++) {
    $aa=$keeps[$ik][5]/$tot1;
    $keeps[$ik][6]=$aa;
    $cum1+=$aa;
    $keeps[$ik][7]=$cum1;

  }

  $ncdf=array_column($keeps,7);  // used in 'First pull'  (pulling from "approximate pdf")
  $retvals=[];
  for ($ido=0;$ido<$ndo;$ido++) {
    $jtrys=0   ;
    for (;;) {    // importance sample  -- keep going until success
      $jtrys++;
      $gog2=mt_rand()/mt_getrandmax();        // pull a "bin" -- the approximate pdf
      $jj=ss\binarySearchMatrix($ncdf, $gog2 )+1    ;
      $jj=max(0,$jj);
      $xmin=$keeps[$jj][0];$xmax=$keeps[$jj][2];$pmax=$keeps[$jj][4];
      $dx=$xmax-$xmin;
      $gog3=mt_rand()/mt_getrandmax();         // find a value in the bin
      $xt=$xmin+($dx*$gog3);
      $xt=intVal(round($xt,0,PHP_ROUND_HALF_UP));
      $p2=binomPdf($nTrials,$xt,$p);
      $gog4=mt_rand()/mt_getrandmax();          // the importance sampling step!
      $pcomp=($p2/$pmax);
      if ($gog4<$pcomp) {                   // passes imporant sampling... so use it
         $retvals[]=$xt;
         break;                        // succes. so go  get next rv
      }
    }   // got 1 rv
  }     // got ndo
  return $retvals;
}


