<?php
namespace wSurvey\Math\rbinomial ;

require ('logGamma.php' );     // must be in same directory as this fiel
require ('sortSearch.php' );
require ('wMean.php' );

use wSurvey\Math\logGamma as lnG ;
use  wSurvey\Math\sortSearch  as ss ;
use  wSurvey\Math\wMean  as wcmp ;


// 22 April 2021
// reverse binomial random variable generator
// A reverse binomial pdf, rbin(pRate|m,N) returns the  probability of a given pRate given m successes in N trials
// The rv generator generates pRates by randomly drawing from this pdf
//
// Setup.
// Assume that ws_rbinomial.php,  sortSearch.php and logGamma.php  are in the xx path
// Include in your .php file
//
// require_once ( 'ws_rbinomial.php' );
// use \wSurvey\math\rbinomial as lRbin ;
//
//
// Thus, if xx=d:/wamp/www/wsurvey/myProgs
// is where your php file is in
// and it uses files in
//   d:/wamp/www/wsurvey/myProgs/lib
// then
//    require_once($xx.'lib/ws_rbinomial.php'');
// points to
// d:/wamp/www/wsurvey/myProgs/lib/ws_rbinomial.php
//    and these two other files must also exist: be in the same directory
// d:/wamp/www/wsurvey/myProgs/lib/logGamma.php and    d:/wamp/www/wsurvey/myProgs/lib/logGamma.php
//
//Usage (assuming the use statement above was specified):
//   a) lRbin/rbinom_rv($nRvs,$m,$N) --
//            return nRvs "pRate" from a reverse binomial with m success and N trials
// Or, if you want to change some of the default parameters used in internal computations....
//  b)   lRbin/rbinom_rv($nRvs,$rbinSpecs)
//  where $rbinSpecs is created by calling
//     rbinSpecs=  lRbin/rbinom_rv_init($m,$N,$arg1,...)
//   where the arg1,... after $N can be left out (defaults will be used)
// The available options are described below
//
// The second alternative is useful if you are generating a number of sets of rvs from the same (same m and N) reverse binomial distribution
//
//
//=====================================
// generate an array (of length nRvs) of  0.0  to 1.0 pRates -- where the probablity of choosing a pRate (for inclusion in the array
//  is based on the reverseBinomial formed given m successes in N trials
//   $nRvs : Number of random variables, pulled from the specified reverse binomial
//  $v1 : either the rbinSpecs (from a call to rbinom_rv_init. Or m (number of successes);
//  $v2 : if v1=m, this must be N (# of events). If v1 is rBinSpecs -- this is ignored.
//  $simple: if 1, then do NOT use discrete distribution -- simple importance sampling. Note that simple improtance sampling is used if method=2 (in _init)
// returns object with these fields:
//         'rvs' : array of rvs selected from the reverse binomial specificaton (specified by v1 and v2)
//         'nTrys'  : number of attempts need to finds these rvs.  If >= $maxtries, count($rvs) usually < $nRvs
//         'nbad' : number of bad attempts -- importance sample failure due to reverse binomial pdf of candidate rate > "discrete distribution"  pdf
//                    If > 0, probably should use a differnt "discrete distribution" (re call rbinom_rv_init)
//
// Notes:
//     If maxTries (default is 10*$nRvs, or you can specify in rbinSpces) is close to nRvs, or the specification is troublesome
//      (a low m/N and a low N), rvs might not have nTrys rows --


function rbinom_rv($nRvs,$v1,$v2=0,$simple=0) {

   $nargs=func_num_args() ;   // if 2, then $v1 must of been creaetd by rbinom_rv_init
   if ($nargs==2 || is_array($v1) ) {

       if (!is_array($v1) || !array_key_exists('ok',$v1) || $v1['ok']!=1) return 'Error: not initialized properly: '.$v1['ok']  ; // an error
       $rbinSpecs=$v1;
       $maxTries=$rbinSpecs['maxTries'] ;
       if (!is_numeric($maxTries) || $maxTries==0 ) $maxTries=10*$nRvs ;
   } else {              // v1=m, v2=N  -- create rbin object with defaults -- such as method=1 (discreteDist method)
      if ($v1>$v2)     return "Error: unable to initialize properly: (m) $v1 lt  (N) f$v2"   ;
      $rbinSpecs=rbinom_rv_init($v1,$v2);
      if (!array_key_exists('ok',$rbinSpecs) || $rbinSpecs['ok']!=1) return "Error: unable to initialize properly with $v1, $v2"   ; // an error
      $maxTries=10*$nRvs;     // default number of maxtries (based on nrvs)

   }
   $m= $rbinSpecs['m']; $N=$rbinSpecs['N'];
   if ($m<0 || $N<0) {
      $allValsI=array_fill(0,$nRvs,0.0);
      return ['rvs'=>$allValsI,'ntrys'=>$nRvs,'nbad'=>0,'message'=>"m=$m, N=$n"];
   }
   if ($m>$N) {
      $allValsI=array_fill(0,$nRvs,1.0);
      return ['rvs'=>$allValsI,'ntrys'=>$nRvs,'nbad'=>0,'message'=>"m=$m, N=$n"];
   }


  $discreteDist=$rbinSpecs['discreteDist'];

  $nbad=0;   // # of unallowed importance sample (pdf of pRate > discreteDistributon  pdf at pRate
  $rbobjUse=$rbinSpecs['rbobj'];  // rbobj   is used to quickly create a binomial pdf value (for fixed m and N)

  $ntry=0;
  $allValsI=[];
  $nsucc=0;$nfails=0;
  $rmethod=$rbinSpecs['method'];


// minimal info -- giveu p
if ($rmethod==0) {
       return ['rvs'=>[],'ntrys'=>0,'nbad'=>0,'message'=>"minimal information initialization: unable to generate rvs "];
}

// simple -- dont use discrte distribution
  if ($simple==1 || $rmethod==2) {               // if method =2, no discreteDist info, so have to use simple
     $topVal=$rbinSpecs['obsRatePdf']*1.1;
     $x1=max(0.0,$rbinSpecs['rangeLowest']/$rbinSpecs['xMinMult']);
     $x2=min(0.999999,$rbinSpecs['rangeHighest']*$rbinSpecs['xMaxMult']);
     $xrange=$x2-$x1;

     while (1) {
        $ntry++;
        $arand=mt_rand()/ mt_getrandmax();
        $xtest=$x1+($arand*$xrange);
        $xtestPdf=rbinomPdf($rbobjUse,$xtest);        // actual pdf of reverse binomial at this rate
        if ($xtestPdf>$topVal) {                     // should not happen ... this is not allowed (imprtance sampling will under choose this point
           $nbad++;
           continue;
       }
        $arand3=mt_rand()/ mt_getrandmax();
        $impTest=$arand3*$topVal;   // the test ....
        if ($impTest<$xtestPdf) { // use it
            $allValsI[]=$xtest ;
            $nsucc++;
            if (count($allValsI)>=$nRvs) break;  // done!
        }  else {
          $nfails++;
        }
        if ($ntry>$maxTries) break  ; // giveup
     }  // infinite loop
     $actualTries=$ntry ;
     return ['rvs'=>$allValsI,'ntrys'=>$ntry,'nbad'=>$nbad,'message'=>"simple importance sampling:   $nRvs : got $nsucc  in $actualTries (fails=$nfails, bad=$nbad)"];
  }   // simple

// ------- cdfs

if ($rmethod==3) {
 
  $cdflist=$rbinSpecs['cdfList'] ;

  for ($ix=0;$ix<$nRvs;$ix++)  {
     $arand=mt_rand()/ mt_getrandmax();
     $rvBin=ss\binarySearchMatrix($cdflist,$arand,2)  ;     // what "bin"  of the cdfList contains this  == add 1 since cdflist is "value and end of bin"
     if ($rvBin>=count($cdflist)-1) {    // last value -- should never happen
         $allValsI[]=$cdflist[count($cdflist)-1][0];   // use highest in range
         $nsucc++;
         continue;
     }
     if ($rvBin<0)  {     //before cdf of first row
         $x0=0; $y0=0;
         $x1=$cdflist[0][0] ; $y1=$cdflist[0][2] ;
     } else {  // inbetwee
         $x0=$cdflist[$rvBin][0];    $y0=$cdflist[$rvBin][2] ;
         $x1=$cdflist[$rvBin+1][0];   $y1=$cdflist[$rvBin+1][2];
     }

     $xrange= $x1-$x0 ;
     $yrange= $y1-$y0 ;
     $dy=$arand-$y0;
     if ($yrange==0)     {  // might happen
        $allValsI[]=$cdflist[$x0][0];   // assign to x0
         $nsucc++;
        continue;
     }
     $dyFrac=$dy/$yrange;
     $xuse=  $x0+($dyFrac*$xrange) ;
     $allValsI[]= $xuse ;   // assign to x0
     $nsucc++;
  }
  return ['rvs'=>$allValsI,'ntrys'=>$nsucc,'nbad'=>0,'message'=>"CDF draws of $nRvs values "];

}

// --- if here, use discreteDist augmented importance sampling
   $nfails=0;$nsucc=0;

   while (1) {
    $ntry++;
     if ($ntry>$maxTries) break  ; // giveup

    $arand=mt_rand()/ mt_getrandmax();
    $afind=$arand*$rbinSpecs['cumArea'];
 
// select, using binary search within cumulative cdf one of the bins (between k1 and k2, to shorten region of binary search)

     $rvBin=ss\binarySearchMatrix($discreteDist,$afind,5) ; // ,$k1,$k2)  ;     // what "bin"  of the discreteDistro to look in (from bins k1 to k2) -- speed up search a bit

     $tt=rbinomial_fromBin($rvBin,$discreteDist,$afind );                   // what point in this bin?

     $rvBinPoint=$tt[0] ; $binPdf=$tt[1];                    // the point, and the bin's probability

     $rbinReal=rbinomPdf($rbobjUse,$rvBinPoint);        // actual pdf of reverse binomial at this rate=rvBinPoint

     if ($rbinReal>$binPdf) {                     // should not happen ... this is not allowed (imprtance sampling will under choose this point
        $nbad++;
        continue;
     }
     $arand3=(mt_rand()/ mt_getrandmax());
     $atest=$arand3*$binPdf;   // the test ....

     if ($atest<$rbinReal) {                    // this is the importance sampling
           $allValsI[]=$rvBinPoint ;
           $nsucc++;
     }   else {
        $nfails++;
     }
     if (count($allValsI)>=$nRvs) break;  // done!

 }  // infinite loop
 $actualTries=$ntry ;


 return ['rvs'=>$allValsI,'ntrys'=>$ntry,'nbad'=>$nbad,'message'=>"discreteDist importance sampling: $nRvs : got $nsucc  in $actualTries (fails=$nfails, bad=$nbad)"];
}

//============================================================
// create a "discrete distribution" , and other specs, given m and N.
//  $rbinSpecs=rbinom_rv_init($m,$N,$minPdfBound,$yjumpMax,$maxRecurse,$nSds,$maxTries);
//     Everything after $N is optional
//        minPdfBound  : (Defaultvalue 0.0001) -- used to set the range of support for pRates.
//                Chooses a pRate (lower & upper bound) s.t.  rbin(pRateL;m,N) < minPdfBound
//                This should be small, especially if m/N <0.1  (
//                A "after range of support" bin is added that uses this as the height, and if this is small the bin area can be quite large (in cdf terms)
//                Note: we do NOT recommend using a value less than the default (1e-7).
//
//       yJumpMax : (defalue value=0.015) -- used to fill in values in  a discrete distribution.
//                  The discrete distribution is used to roughly approximate the target continuous reverse binomial pdf.
//                  This discrete distribution is used to efficiently draw rvs (using importance sampling) from the desired reverse binomial.
//                  Points (pRate(i) are added to this discrete s.t.  rbin(pRate_i) - rbin(pRate(i-1)) < yJumpMax ---
//                   where i=0..N are  N points (between bounds of 0.0 and 1.0) that are defined for this discrete distirbution.
//                   Note that for very tight distributions, this may not be achieved -- see the maxRecurse parameter
//        Smaller values of these parameters lead to more points in the discrete distribution. This may lead to fewer attempts used to
//        find the nRvs values. However, each attempt might take longer.
//
//     $maxRecurse : (default=6) how many times to go through the discrete distribution.  Larger values make it more likely to achieve the yJumpMax
//                 target, but will also require longer initialization time.
//     $nSds : (default =2) values to always include in the discrete distribution  obsRate + (and -) 1..nSds standard deviations
//             Using the binomial sd (standard deviation) of m at obsRate and N:  compute and use pRates calculated
//                    as m +/- j* sd ; where (j=1...nSds);
//             More values might help select from tight ( large m and N) rbin distributions. Too large will add unuseful proprocessing time
//             Note: at pRate =0.0 or 1.0; since sd is NOT defined. Hence, nSds is not used (and  1/($N-sqrt($N)) is used to start bounds search).
//
//     $maxTries  : (default= $nRvs*10). If 0, default is used).  Maximum number of attempts to find the desired number of rvs (in
//                  rbinom_rv). If the discrete distribution poorly matches the underlying reverse binomial, importance sampling can fail frequently
//                  This limits the number of attempts. Once this limit is hit, the array of rvs is returned (even if it is length is < $nRvs)
//
//   $minBins : (default=0). If not 0, the minimum numbef of bins -- rows in disceteDist.
//              If the normal algorithim does not yield this many bins, yJumpMax is decrease and another recursion is done.
//
//   $useRates : (default=[]) : array of rates which should be included in discreteDist. If empty, no extra rates adde.
//                 A csv is converted into an array
//
//   $method : default=1 ; 0=test (return rbinomPdf_init info), 1=discreteDist assisted importance sample (the default), 2:raw important sample, 3:cdf lookup
//
//   xMinMult and xMaxMult : (default 1 and 1) not supported in arguments. Supported in altenative calling method.
//                           After finding min and max (that satisfy minPdfBound, multiply the min by xMinMult and the max by xMaxMult
//                                 (constraining to fall within 0.0 and 1.0)
//                           Return these in rbinSpecs (rangeLowest and rangeHighest)
//                           xMinMult and xMaxMult are only used for method=3 (cdf lookup)
//
// returns associative array  that is used by  rbinom_rv.
// One of the fields is 'discreteDist' :an array of "points of support" for the discrete distribution used to pull random variables from
// Each row is a "bin" (one could think of it as a bin in a histogram) -- consisiting of an array with 6 values:
//  [0] - pvalue -- candidate value
//  [1] - rbinomial pdf value (at pvalue)
//  [2] - width of bin (pvalue to pvalue of  next bin)
//  [3] - height of bin  (max of [1] of this and [1] of next bin. For last bin, just use this bins [1]
//  [4] - area of bin (width [2] * [3]
//  [5] - normalized cumulative area  of [4] -- this is used in binary searches (for a value between 0 and 1)
//
// brinSpecs['cdfList'] can also be added -- an 'evenly' (and closely) space matrix of values from the range of support of the rbinomial (
//  [0] : rate
//  [1] : nomralized [df
//  [2] : cumulative cdf
//  [3] : original (non-normalized) pdf
//
// april 2021. Alt  caling (recommended)
//     rbinom_rv_init($m,$N,$options) -- options is assciative array with the arguments (case insensitive)
// if options is specified, the other arguments  (arg 4,...) are ogmpred.
// If an option is NOT specified in options, defautls are used


function rbinom_rv_init($m,$N,$minPdfBound=1e-7,$yjumpMax=0.015,$maxRecurse=30,$nSds=3,$maxTries=0,$minBins=0,$useRates=[],$method=3) {

  $obsRate=$m/$N;

  if (is_array($minPdfBound)) {  // alternative calling method
     $rbinSpecs=['minPdfBound'=>1e-7,'yjumpMax'=>0.015,'maxRecurse'=>10,'nSds'=>3,
              'maxTries'=>1, 'minBins'=>0,'useRates'=>[],
              'method'=>1,'xMinMult'=>1.0,'xMaxMult'=>1.0 ] ;  // the defaults

     $daopts=$minPdfBound;
     $zlookups=[];
     foreach ($rbinSpecs as $aopt=>$goo) {  // use this to support case insensitve lookups
         $lopt=strtoupper($aopt);
         $zlookups[$lopt]=$aopt;
     }
     foreach ($daopts as $aopt=>$optval) {   // find known options, save to rbinSpecs
        $topt=strtoupper($aopt);
        if (!array_key_exists($topt,$zlookups)) continue ; // unknown option
        $saveOpt=$zlookups[$topt];
        $rbinSpecs[$saveOpt]=$optval;
     }

  }  else {   // use the provided arguments, or their default values

     $rbinSpecs=['minPdfBound'=>$minPdfBound,'yjumpMax'=>$yjumpMax,'maxRecurse'=>$maxRecurse,'nSds'=>$nSds,
              'maxTries'=>$maxTries,'minBins'=>$minBins,'useRates'=>$useRates,
               'method'=>1,'xMinMult'=>1.0,'xMaxMult'=>1.0  ] ;  // arguments, or default values of argumetns
  }


  $rbinSpecs['m']=$m;  $rbinSpecs['N']=$N;   $rbinSpecs['obsRate']=null;
  $rbinSpecs['ok']=0;  $rbinSpecs['rbobj']=[];   $rbinSpecs['discreteDist']=[];

  $useRates=$rbinSpecs['useRates'];    // convert csv to array
  if (is_string($useRates)) $useRates=str_getcsv($useRates);
  if (!is_array($useRates)) $useRates=[];
  $rbinSpecs['useRates']=$useRates;

  $method=$rbinSpecs['method'];   // make sure  method is int between 0 and 3
  $method= (!is_numeric($method)) ? 1 : intval($method);
  if ($method<0 || $method>3) $method=1;
  $rbinSpecs['method']=$method;

  if ($m<0 || $m>$N || $N<=0) return $rbinSpecs ; // degenerate cases

  $rbinSpecs['obsRate']=$m/$N ;
  $rbinSpecs['rbobj']=rbinomPdf_init($N,$m);         // this is use by rbinomPdf
  $rbinSpecs['obsRatePdf']=rbinomPdf($rbinSpecs['rbobj'],$rbinSpecs['obsRate']);

// minimal info method ....
  if ($method==0) {        // just return rbobj (and obsrate); and info for useRates
     sort($useRates,SORT_NUMERIC);
     $gotsUse=[];
     foreach ($useRates as $jj=>$arate) {
         if (!is_numeric($arate) || $arate<0.0 || $arate>1.0) continue ;
         $aratePdf=rbinomPdf($rbinSpecs['rbobj'],$arate);
         $gotsUse[]=[$arate,$aratePdf];            // this is the mode
     }
     $rbinSpecs['discreteDist']=$gotsUse ;
     $rbinSpecs['ok']=1;
     return $rbinSpecs ;                   // minimal info
  }

  $ggUse=rbinomial_findBounds($rbinSpecs);  // initialze the discreteDist with obsRate and a series rate (above and below) that find the range of support

  if (count($useRates)>0)  {   // add in "special rates
     sort($useRates,SORT_NUMERIC);
     $gotsUse=[];

     foreach ($useRates as $jj=>$arate) {
         if (!is_numeric($arate) || $arate<0.0 || $arate>1.0) continue ;
         $aratePdf=rbinomPdf($rbinSpecs['rbobj'],$arate);
         $gotsUse[]=[$arate,$aratePdf];            // this is the mode
     }
     if (count($gotsUse)>0) $ggUse=ss\merge_ordered($ggUse,$gotsUse) ;      // quick merge arrays, resulting array ordered on [0] element
  }

// simple importance sampling method .....

  if ($method==2)  {       // no discreteDist important sample. Just uses x range info
  //$rbinSpecs['obsRate']
    $rbinSpecs['rangeLowest']=$ggUse[0][0];
    $rbinSpecs['rangeHighest']=$ggUse[count($ggUse)-1][0];
    $rbinSpecs['ok']=1;
    $rbinSpecs['discreteDist']=$ggUse  ; // ifnormational version (can't be used for augmented importance sampline
    return $rbinSpecs;         // note: xMinMult and xMaxMult rbinom_rv
  }             // end of simple importance sampling method .....

// cdf method .....
  if ($method==3)  {       // no discreteDist important sample. Just return xrange
    $rbinSpecs['rangeLowest']=$ggUse[0][0];    // these may be adjusted in makeceflist
    $rbinSpecs['rangeHighest']=$ggUse[count($ggUse)-1][0];
    $rbinbSpecs['obsRate']=$obsRate;
    $cdfList= rbinomial_makeCdfList($rbinSpecs);
    $rbinSpecs['cdfLowest']=$cdfList[0][0];                // fix it up
    $rbinSpecs['cdfHighest']=$cdfList[count($cdfList)-1][0];
    $rbinSpecs['cdfList']=$cdfList;
    $rbinSpecs['ok']=1;
    return $rbinSpecs;
  }

// ==== if here, method=1 (discretDist)

  $ggNew=rbinomial_fill($ggUse,$rbinSpecs)  ; // fill in with pRates (to reduce jumps in pdf value)

//  $aa1=do_dump_string($ggNew );

  $ggNew=rbinomial_cleanup($ggNew,$rbinSpecs)  ; // add start and end rates, remove small steps

//  $aa=do_dump_string($ggNew );
//  wsretReturn($aa1.'<hr>'.$aa,['status'=>1,'plotInfo'=>[]],1);

  $ggNew=rbinomial_addBinArea($ggNew,$rbinSpecs)  ;  // add in some cdf etc info


  $rbinSpecs['method']=$method;
  $rbinSpecs['discreteDist']=$ggNew ;
  $rbinSpecs['ok']=1;   // signal success

  return $rbinSpecs ;

}


// ==============
// cresat an array of rbinomial rate,pdf,cdf rows; equally spaced between xmin and xmax (adjusted by xminMult and xmaxMult
// spacing uses minObs (if not 0). Otherwise uses yjumpMax (as a fraction of the obsRate
// modified rbinSpecs as call by ref
// each row: aRate, normalizedPdf, cdf (cdf = 1 in last row), originalPef
function  rbinomial_makeCdfList($rbinSpecs) {
  
  $foores=['status'=>1,'plotInfo'=>[]];
    $robj=$rbinSpecs['rbobj'];

    $obsRate=$rbinSpecs['obsRate'];
    $rangeLowest=max(0,$rbinSpecs['rangeLowest']*$rbinSpecs['xMinMult']);
    $rangeHighest=min(0.999999,$rbinSpecs['rangeHighest']*$rbinSpecs['xMaxMult']);
    $xRange=$rangeHighest-$rangeLowest;
    $nX=$rbinSpecs['minBins'];
    $yjump=$rbinSpecs['yjumpMax'];

    if ($obsRate==0)   {          // special case
      $adenom=1/sqrt($rbinSpecs['N'] );
      $pdf0= rbinomPdf($robj,0.0);  // should be 1.0
      $amult=0.05;
      $tryRate=$adenom;

      while (1)    {          // search for small enough xinc..
          $pdf1= rbinomPdf($robj,$tryRate);
          $adiff=abs($pdf0-$pdf1);
          if ($adiff< $yjump) {
             $xinc=$tryRate;
              break;
          }
           $tryRate=$amult*$tryRate;
       }

    } else if ($obsRate==1.0)   {  // other special case
          $adenom=1/sqrt($rbinSpecs['N'] );
         $pdf0= rbinomPdf($robj,1.0);  // should be 1.0
         $amult=0.05;
         $tryRate=$adenom;

        while (1)    {          // search for small enough xinc..
         $tryRate1=1-$tryRate;
          $pdf1= rbinomPdf($robj,1-$tryRate);
          $adiff=abs($pdf0-$pdf1);
          if ($adiff< $yjump) {
             $xinc=$tryRate;
              break;
          }
          $tryRate=$amult*$tryRate;

       }

    } else {
       if ($nX>0)  {
          $xinc=$xRange/$nX;
       } else {
         $xinc=$obsRate*$yjump ;
      }
    }

 
    $nsteps=$xRange/$xinc;
    $crates=[$obsRate];
    $ado1=$obsRate  ;
    $ado2=$obsRate ;


    while (1)  {
       $idid=0;
       $ado1-=$xinc;
       if ($ado1>=$rangeLowest) {
          $idid++;
          array_unshift($crates,$ado1);
       }
       $ado2+=$xinc;
       if ($ado2<=$rangeHighest)  {
          $crates[]=$ado2;
          $idid++;
       }
       if ($idid==0 || count($crates)>10000) break ;   // at the low and high bounds. Or very big
    }       // got an evenly seperate list of ranges between xlow and xhigh, that always is centered on obsrate


    $gotsUse=[];
    for ($nn=0;$nn<count($crates);$nn++) {        // compute rbin probablity at each of these rates
       $arate=$crates[$nn];
       $aratePdf=rbinomPdf($robj,$arate);
       $gotsUse[]=[$arate,$aratePdf,0,$aratePdf];
    }
    $acum=0;
    

    for ($nn=0;$nn<count($gotsUse);$nn++) $acum+=$gotsUse[$nn][1] ;
    for ($nn=0;$nn<count($gotsUse);$nn++)  $gotsUse[$nn][1]/=$acum ;  // normalize to pdfs add to 1
    $acum=0;
    for ($nn=0;$nn<count($gotsUse);$nn++) {     // compute running cdf at each of these rates
       $acum+= $gotsUse[$nn][1];
       $gotsUse[$nn][2]=$acum;
    }


    return $gotsUse ;
}


//=====================      
// find upper and lower bounds (retaining ssome stuff inbetwee
function rbinomial_findBounds($rbinSpecs ) {
   $ggUse=[];
   $rbobj=$rbinSpecs['rbobj'];   $minPdfBound=$rbinSpecs['minPdfBound'];

  $ggLow=rbinomial_findBound(0,$rbobj,$minPdfBound);
  $ggUse=array_reverse(array_slice($ggLow,1));      // renove the "current rate" (its retainged in gghHigh
  $ggHigh=rbinomial_findBound(1,$rbobj,$minPdfBound);
  foreach ($ggHigh as $ii=>$gval) $ggUse[]=$gval;    // merge high bounds list after low bounds  list

  return $ggUse ;
}

//=================
// find lower, or upper, bound  of rbin -- the "p" value at which the reverse binomial (given N and m) is small enough to ignore
//  adire: -1 : lower bound, 1=upper bound
//  $rbobj : object created by rbinomial_init (contains N,m, and other values to use in a binomial pdf calc)
//  $minPdfBound : minimum value of pdf. Keep lookihng (increasing or decreasing the bounds) for a p value whose rbinomial pdf value is < minPdfBound
// note: assume that 0<m<N

function rbinomial_findBound($adire,$rbobj,$minPdfBound) {
   $gots=[];
    $N = $rbobj['N'];          // # events
    $m = $rbobj['m'];          // prob of success

   $pTest= $m/ $N  ;

   $pTestPdf=rbinomPdf($rbobj,$pTest);

   $gots[0]=[$pTest,$pTestPdf];            // this is the mode

   if ($adire==0) {                        // lower bound
     if ($m==0) {            //degenerate case  -- cant get lower
        return $gots;
    }
    $amult1=0.86 ; $amult2=0.66;
   } else {       // upper bound
      if ($m==$N) {            //degenerate case, cant get higher
        return $gots;
      }
      $amult1=1.14 ; $amult2=1.33;
      if ($m==0) $pTest=1/($N-sqrt($N))  ; // avoid trap at success rate=0
   }

   $pTest=$amult1*$pTest;                 // always have at least one step
   $pTestPdf=rbinomPdf($rbobj,$pTest) ;
   $gots[1]=[$pTest,$pTestPdf];
   $indGot=1;

   for ($ii=0;$ii<100;$ii++) {    //  keep looking below --  until most recent pdfMeasure is < $minPdfBound.
       if ($pTestPdf<$minPdfBound  )  break;  // small enough!
       $indGot++ ;
       $pTest=$amult2*$pTest;
       $pTestPdf=rbinomPdf($rbobj,$pTest) ;
       $gots[$indGot]=[$pTest,$pTestPdf];
   }


  return $gots ;

}

//=================
// recursively insert rows  in array of ordered rbin pairs (pval,pvalPdf). Add pairs s.t. there aren't large jumpts in the y value (the pdf value)
// eg:  abs(pvalPdf[n]-pvalPdf[n-1]) < yjump)
// returns ordred (in pval) array

function rbinomial_fill($gvals,&$rbinSpecs,$ict=0,$yjumpMaxUse=0)   {
  $gadds=[];
  $rbobj=$rbinSpecs['rbobj'];
  $yjumpMax= ($yjumpMaxUse<=0.0) ? $rbinSpecs['yjumpMax'] : $yjumpMaxUse ;     // use default in rbinspecs, but allow for override
  $maxRecurse=$rbinSpecs['maxRecurse'];
  $nsds=$rbinSpecs['nSds'];
  $minBins=$rbinSpecs['minBins'];
  if ($ict==0)   {      // SPECIAL CASE: on first call, add some reasonable values (using sd)
      $rbinSpecs['irecurse']=0;   // add in a few values using sd of binimoa
      $rbinSpecs['yjumpMaxUse']=$yjumpMax;   // might be changed if minBins>0

      if ($nsds>0) {
        $pObs=$rbinSpecs['obsRate'];
        $N=$rbinSpecs['N']; $m=$rbinSpecs['m'];
        $aSd=sqrt($N * $pObs * (1-$pObs));       // this does nothing if m=0 or m=N
        if ($aSd==0) $aSd=max(1,log($N));
        $asd2=$aSd*0.5 ; $topsd=$aSd*$nsds;
        $useSd=-($aSd *$nsds) ;
        while (1) {
          if ($useSd==0) {
             $useSd+=$asd2;
             continue;
          }
          $pRate=($m+$useSd)/$N;
          if ($pRate<=0.0) {
             $useSd+=$asd2;
             continue;
          }
          if ($pRate>1.0) break;
          $pratePdf=rbinomPdf($rbobj,$pRate);
          $gadds[]=[$pRate,$pratePdf];
          $useSd+=$asd2;
          if ($useSd>$topsd) break;
        }
//        $gvals=ss\merge_ordered($gvals,$gadds) ;
       $gvals=\wSurvey\Math\sortSearch\merge_ordered($gvals,$gadds) ;
        $gadds=[];    // reset
      }

  }
// look at adjacent pairs in gvals. If they are "far apaert" (differnce in pdfs), add a point between them
  $ptestPdf=$gvals[0][1];

  for ($ii=1;$ii<count($gvals); $ii++) {
     $priorPdf=$ptestPdf;
     $ptestPdf=$gvals[$ii][1];
     $adiff=abs($ptestPdf-$priorPdf);
     if ($adiff<$yjumpMax) continue  ;  // not too far apart

     $priorP=$gvals[$ii-1][0];
     $pnow=$gvals[$ii][0];
     $pnew=$priorP + (($pnow-$priorP)/2) ;
     $pnewPdf=rbinomPdf($rbobj,$pnew);
     $gadds[]=[$pnew,$pnewPdf];
  }

  if (count($gadds)==0)  {       // nothing needed to be added -- done? Check minbins!
     if (count($gvals)>=$minBins || $ict==$maxRecurse)    return $gvals ; // nothing added
     $yjumpMaxUse=0.66*$yjumpMax  ;      // need more bins .. try with tighter bounds
     $ict++;
     $rbinSpecs['yjumpMaxUse']=$yjumpMaxUse;
     $rbinSpecs['irecurse']=$ict;
     $ggNewB=rbinomial_fill($gvals,$rbinSpecs,$ict,$yjumpMaxUse) ;
     return $ggNewB ;  // got enough bins
  }

// if here, have some additions. Insert them maintaining order

//  $ggNew=ss\merge_ordered($gvals,$gadds) ; // merge -- ggNew is sorted by [0] (pval)
  $ggNew=\wSurvey\Math\sortSearch\merge_ordered($gvals,$gadds) ; // merge -- ggNew is sorted by [0] (pval)

  $maxDiff=0;  // compute max differnce -- if too big, recurse to addin some more
  for ($ig1=1;$ig1<count($ggNew);$ig1++)  {
     $adiff=abs($ggNew[$ig1][1]-$ggNew[$ig1-1][1]);
     $maxDiff=max($adiff,$maxDiff);
  }
  if ($ict<$maxRecurse) {           // not too many recursions
    if ($maxDiff>$yjumpMax  || count($ggNew)<$minBins) {       // recurse to fill in some more -- at some point stop, ie if a lot of area on one point
     $ict++ ;
     $rbinSpecs['irecurse']=$ict;
     $ggNew2=rbinomial_fill($ggNew,$rbinSpecs,$ict) ;      // a bit inefficient, but good enough for now (might waste time checking maxdiff
     return $ggNew2 ;
    }
  }


  return $ggNew ;  // good enough
}


//=================
// add to ends of list of rates, remove redundant rates, remove too close rates
 function rbinomial_cleanup($ggNew,&$rbinSpecs)  { // add start and end rates, remove small steps
  $ggNew1=[];         // remove duplicates
  $ilast=count($ggNew)-1;

  $yjump= ($ilast<20) ? 0 :  $rbinSpecs['yjumpMax']/10 ;
 
  $oldval=$ggNew[0][0];
  foreach ($ggNew as $igg=>$agg) {      // remove duplicates
      $arate=$agg[0];
      if ($arate!=$oldval) $ggNew1[]=$agg;
      $oldval=$arate ;
  }

  $ilast=count($ggNew1)-1;

  if ($rbinSpecs['minBins']<=0) {   // if minbins specified, don't clear out interiors
    $ggNew2=[];  // remove too close
    $ggNew2[0]=$ggNew1[0] ; $aval=$ggNew1[0][1]; $lastSave=0;
    $ij=1;
    while (1) {
      if ($ij>$ilast) break;   // at end of list
      $aval2=$ggNew1[$ij][1];
      $adiff=abs($aval2-$aval);
      if ($adiff>$yjump)  {    // ydiff now too big. Save most recent and reset stuff
          if ( ($ij-$lastSave)==1 ) {   // nothing to remove
             $ggNew2[]=$ggNew1[$ij];
             $lastSave=$ij;
             $aval=$aval2;
             $ij++;
             continue;
          }
          $ggNew2[]=$ggNew1[$ij-1];     // the last point with yjump stiill > diff
          $aval=$ggNew1[$ij-1][1];
          $lastSave=$ij-1;  // do NOT update ij (so look at ij again, comparing to prior, not to "most recent saved"
          continue ;
       }   // else got to next (do not copy ij just yet, perhaps never.

       $ij++ ;         // look at next one
    }       // while
    if ($lastSave!=$ilast) $ggNew2[]=$ggNew1[$ilast];
  } else {                                           // no interior clearout
    $ggNew2=$ggNew1;
  }

  $ilast=count($ggNew2)-1;
  $rbinSpecs['rangeLowest']=$ggNew2[0][0];    // update these "non-terminal" ranges
  $rbinSpecs['rangeHighest']=$ggNew2[$ilast][0];
  $rbinSpecs['nRange']=count($ggNew2);   // the "supported" x range (outside of this covered by prepended and appended bins.

// add widths
  $rbobj=$rbinSpecs['rbobj'] ;

// add the 0.0 and 1.0 actual bounds
   $ilast=count($ggNew2)-1;

// first bin: 0 to start of x support. Skip if already specified
  if ($rbinSpecs['rangeLowest']>0.0) {
     $rate1=$ggNew2[0][0];
    if ($rbinSpecs['m']==0) {
      $add1=[0.0,1.0];      // add to start. Use 1 as the pdf, since m=0
    } else {
      $add1=[0.0,0.0];      // add to start. Use 0 as the pdf, since it is never used
    }
 }

// last bin : last of x suport to 1.0
  if ($rbinSpecs['rangeLowest']<1.0) {
    $ilast=count($ggNew)-1;        // the last "range of support bin"
    $arate=$ggNew[$ilast][0]*  1.000001  ;   // use   a rate just past the "last supported x"
    $pnewPdf=rbinomPdf($rbobj,$arate);         // its rate
    $add2=[$arate,$pnewPdf,1.0-$arate,$pnewPdf];      // rate and   pdf and binWidth  and binHeight
  }


  array_unshift($ggNew2,$add1);  // prepend
  $ggNew2[]=$add2;  // append


// add widths (but not to   1.0 rows which has a width already)
  for ($ij=0;$ij<count($ggNew2)-1 ; $ij++) {  // last one has special rate set above (1-...)
     $ggNew2[$ij][2]=$ggNew2[$ij+1][0]-$ggNew2[$ij][0];
     $ggNew2[$ij][3]=max($ggNew2[$ij+1][1],$ggNew2[$ij][1]);

  }


  return $ggNew2 ;
}

//==============
// modifies the discreteDist -- adds "bin areas"
//  Each row of this array can be thought of a bin of a discrete distibution (a histogram shaped pdf that approximates the true reverse binomial pdf)
//  Each row has:
//  [0] - pvalue
//  [1] - rbinomial pdf value (at pvalue)
//  [2] - width of bin (from left size at pvalue to right side of the pvalue of the next bin)
//  [3] - height of bin  (max of this, and the next bin's, rbinomial pdf. Saved for informational purposes -- is use to compute impScale
// and added here...
//  [4] - area of bin (width x height)
//  [5] - running sum of bome area -- this is used in binary searches (for a value between 0 and 1)
//
//
// Also sets some rbinSpecs parameters (call as ref)
//
// Note: [1] is the "pdf" at the pvalue (the left end of hte bin
//       [4] is the "area" of the bin -- kind of a "discrete" pdf"
//       [5]  is the "cdf" of the "right edge" of the bin. So -- it is actually the cdf at the pvalue of the NEXT bin.
//            For the last element, the next bin's pvalue is 1.0.

function rbinomial_addBinArea($ggNew,&$rbinSpecs) {   // rbinSpecs changed in place


// now start to add  width,height, and area info to bins in the  range of support
  $iModeAt=0;
  $modeY=0;

  $yMin=11; $yMax=-1; $xMin=11; $xMax=-11 ; $yMaxI=0;

// do stuff in the "range of support" (don't do the stuff  prepended and appended
  for ($ig1=0;$ig1<count($ggNew);$ig1++)  {    // add bin width/height/area .
    $arate=$ggNew[$ig1][0];
    $apdf=$ggNew[$ig1][1] ;  $yMin=min($yMin,$apdf); $yMax=max($yMax,$apdf);
    if ($apdf==$rbinSpecs['rangeHighest'])$yMaxI=$ig1;
    if ($ggNew[$ig1][0]==$rbinSpecs['obsRate']) {   // the observed rate? Save some extra info
       $rbinSpecs['obsRateIndex']=$ig1;
       $rbinSpecs['obsRatePdf']=$apdf ;
    }
    $binArea=$ggNew[$ig1][2] * $ggNew[$ig1][3] ;  // area in this bin
    $ggNew[$ig1][4]=$binArea;

  }
  $rbinSpecs['minRangePdf']=$yMin;
  $rbinSpecs['maxRangePdf']=$yMax;
  $rbinSpecs['iMaxRangePdf']=$yMaxI ;


// compute normzlied bin area (a "bin pdf") and bin area cdfs -- include prepended and appended bins
 $cumArea=0;
  for ($ig1=0;$ig1<count($ggNew);$ig1++)    {
    $binArea =$ggNew[$ig1][4];               // bin area -- used in importance phase
    $cumArea+=$binArea ;
    $ggNew[$ig1][5]=$cumArea ;             // sum of bin area - used in binary search phase
   }
 $rbinSpecs['cumArea']=$cumArea;


 $rbinSpecs['discreteDistComment']='0:rate,1:ratePdf,2:binWidth,3:binHeight,4:binArea,5:binAreaCum';
 return $ggNew ;


}

//================
// randomly choose a point in a bin (from left side to right side of bin
// [0] element is pvalue=0, but cdf from pvalue to [1] element
//  [0] - pvalue
//  [1] - rbinomial pdf value (at pvalue)
//  [2] - width of bin (from left size at pvalue to right side of the pvalue of the next bin)
//  [3] - height of bin  (max of this, and the next bin's, rbinomial pdf. Saved for informational purposes -- is use to compute impScale
//  [4] - area of bin (width x height)
//  [5] - running sum of bome area -- this is used in binary searches (for a value between 0 and 1)
//

function  rbinomial_fromBin($rvBin,$ggNew,$afind ) {
   if ($rvBin<0)  {                            // from 0 to first in ggNew
      $x0=0;  $xrange=$ggNew[0][2];
      $y0=0;      $y1=$ggNew[0][5];
      $aheight=$ggNew[0][2] ;
  } else {
     $x0=$ggNew[$rvBin+1][0]; $xrange=$ggNew[$rvBin+1][2];
     $y0=$ggNew[$rvBin][5];    $y1=$ggNew[$rvBin+1][5];
     $aheight=$ggNew[$rvBin+1][3] ;
  }
  $yrange=abs($y1-$y0);
  $dy=abs($afind-$y0);     // in the "cumulative area" array
  if ($yrange>0) {
    $dx=$dy/$yrange;
    $xfig=$x0+($dx*$xrange);
  } else {
      $xfig=$x0+($xrange/2); //  a hack
  }
//  $say1=" $afind $y0 $y1  :   :  dy=$dy  / yrange=$yrange = dx , $dx * $xrange  + $x0 ==> $xfig ::::  ";
 return [$xfig,$aheight];

}

//========================



//=================
// returns weighted statistics associated with a "discrete reverse binomial distribution" as created by =rbinom_rv_init
//rbinSpecs['discreteDist'][0..nrows]
//  [0] - rate  (a candidate "success rate"
//  [1] - rbinomial pdf value (at pvalue)
//  [2] - width of bin (pvalue to pvalue of  next bin)
//  [3] - height of bin  (max of rbinomial pdf value and the rbinomial pdf value of the next bin
//  [4] - area of bin (width x height)
//  [5] - cumulative area -- this is used in binary searches
// The weights are normalized values of [4] -- the "bin area" . If [4] does not exist (eg; a cdflist is being used), use [1] (the pdf at a rate)
// stats returned (as fields of an associative array)
//  nrows: # or rows
//  rateMin, rateMax: min and max of rates ) -- [0] from discreteDist
//  binMin binMax: min and max of "bin weight", after normalization. (using [4] or [1] from discretDist)
//  binCdfs : array of "cdf" values of the normalized weights. Value of last row is 1.0
//  totBinWeight: sum of weights BEFORE normalization. This is used to normalize!
//  rateMean  rateSd : mean and sd of rates
//  customMean customSd: mean and sd of custom function. If no custom fuction, equals false
//  rates, weights, cdfs, customs: arrays  of  the rates, normalized weights, cdfs (of normalized rates), and the custom values
//
// Custom values are created using the funcname, a candidate rate, and the "non-changing"  parameters from funcargs
// Thus, funcname is used as:
//  customs[i]=funcname(rates[i],funcargs[0],funcargs[1],...);
// Note that the weight is NOT passed to funcname.
//

function rbinom_rv_stats($arbinSpecs,$funcname=false,$funcargs=[]) {

   $rrates=[];$customs=[];$binWeights=[];
   $qcustom=false;

   $stats=['nrows'=>false,'xMean'=>false,'xSd'=>false,'customMean'=>false,'customSd'=>false] ;

   $discreteDist=$arbinSpecs['discreteDist'];
   $ndisc=count($discreteDist);
   $stats['nrows']=$ndisc;
   $stats['xmin']=$discreteDist[0][0];
   $stats['xmax']=$discreteDist[$ndisc-1][0];
   if ($funcname!==false && is_callable($funcname)) {
       $qcustom=true;
       array_unshift($funcargs,0)    ; // a candidate rate  is always first arg to funcname
   }


   $totBinWeight=0;   $binMax=-11111111111;   $binMin=11111111111 ;    $oldcdf=0;
   $cdfs=[];
   for ($ii=0;$ii< $ndisc;$ii++) {
       $rrates[$ii]=$discreteDist[$ii][0] ;
       $abweight= (array_key_exists(4,$discreteDist[$ii])) ? $discreteDist[$ii][4] : $discreteDist[$ii][1]  ;
       $binWeights[$ii]=$abweight ;
       $totBinWeight+=$binWeights[$ii]   ;   // the normalization
       if ($qcustom) {
           $funcargs[0]=$rrates[$ii] ;
           $customs[$ii]=call_user_func_array($funcname,$funcargs);
       }
   }
   $stats['totBinWeight']=$totBinWeight;

   $stats['rates']=$rrates;
   if ($qcustom) $stats['customs']=$customs;

   for ($ii=0;$ii< $ndisc;$ii++) {
       $abweight=$binWeights[$ii]/$totBinWeight ;
       $binMax=max($binMax,$abweight) ;
       $binMin=min($binMin,$abweight) ;
       $binWeights[$ii]= $abweight   ;   // the pdf normalization
       $oldcdf+=$abweight;
       $cdfs[$ii]=$oldcdf;
   }
   $stats['weights']=$binWeights;
   $stats['binMin']=$binMin;
   $stats['binMax']=$binMax ;
   $stats['binCdfs']=$cdfs;

   $s1=wcmp\wcmp_mean_var($rrates,$binWeights);
   $stats['rateMean']=$s1[0]; $stats['rateSd']=sqrt($s1[1]);
   if ($qcustom)  {
       $cust1=wcmp\wcmp_mean_var($customs,$binWeights);
        $stats['customMean']=$cust1[0]; $stats['customSd']=sqrt($cust1[1]);
   }


   return $stats;
}

//--------------------
// probablity function


//================
// reverse binomial init: object contains m and N, returns pdf given p

function rbinomPdf_init($N,$m) {
 $errs=[] ;
 $rets=['N'=>$N,'m'=>$m,'lngN'=>null,'lngm'=>null,'lngNm'=>null,'term1'=>null,'errs'=>$errs ];
 if (is_numeric($N)) {
     $N=intVal($N);
     if ($N<0) {
         $errs[]='N (events) is lt 0 ';
     } else {
         $rets['lngN']=lnG\logGamma($N+1);
     }
  } else {
     $errs[]='N is non-numeric' ;
  }

   if (is_numeric($m)) {
     $m=intVal($m);
     if ($m<0) {
         $errs[]='m (successes) is lt 0 ';
     } else {
         $rets['lngm']=lnG\logGamma($m+1);
     }
  } else {
     $errs[]='m is non-numeric' ;
  }

  if (count($errs)==0){
      $rets['lngNm']=lnG\logGamma(1+$N-$m);
      $rets['term1']= $rets['lngN']-  ($rets['lngm'] +    $rets['lngNm']);
  }
  $rets['errs']=$errs;
  return $rets;
}
// ----- returns pdf value of a  reverse binomial, given an object with N and m, and a value of p.
//  xx is created by rbinomPdf_init(N,m) -- it contains intermediate parameters that do not change with p

 function  rbinomPdf($xx,$p) {             // no 3 arg version (you can use binomPdf for that)

    $N = $xx['N'];          // # events
    $m = $xx['m'];          // prob of success
    
// some edge dases
    if ($p<=0.0) {
       if ($m==0) return 1.0;
       return 0.0;
    }
    if ($p>=1.0) {
       if ($m==$N) return 1.0;
       return 0.0;
    }
// non edge cases
//    $lngN=$xx['lngN'];     // useful term. Note that $m=# of successes
//    $lngm=$xx['lngm'];     // useful term. Note that $m=# of successes
//    $lngNm=$xx['lngNm'] ;
//   $logBinom =  $lngN- $lngm - $lngNm ;
    $logBinom=$xx['term1'];
    $p=floatVal($p);
    $lnp=log($p);         // no object, so calculate these
    $lnp_1=log(1-$p);
// now do the pdf calculation
   $logP = $logBinom + $m * $lnp + (($N - $m)* $lnp_1);

   $vv=exp($logP);
   return $vv ;
}

//====================
// return 1st and second derivative of the rverse binomial, at each point in a array or matrix
//  reverse binomial= pdf =f(rate;N,m) -- N and m are "fixed", what changes is the the underlying "rate"
// Uses standard first and 2nd deriv of a binomial. For discussion see:
// https://math.stackexchange.com/questions/2216306/derivative-of-binomial-probability
// https://math.stackexchange.com/questions/651976/second-derivative-of-binomial-distribution
//
// varray: matrix: each row is an array [0]=rate, [1]= binomial prob at this rate (
//         Or, an array of rates.
//           If an array (and not a matrix), the binomial probablity - a N,m, and the rate in this row -- is calculated.
//           So providing a matrix is a short cut -- but should be used only if you already have the probablity  (at each rate) calculated
// m, N : m successes in N attempts -- not use of prob is provided in varray
// iProb: optional  (default=1). col in matrix containing probablities to compute derivative for. Ignored if varray is not a matrix
// iProb: optional  (default=0). col in matrix containing rates to compute derivative for. Ignored if varray is not a matrix
//
// returns array: [rate,1st deriv, 2nd deriv) (derivs are wrt p, at each rate and using the same N and m for all rows)
//
// first and second  derivative of binmoial f(p|N,m). Fixed n and m, array of p values
// f(p)= bin(p|m,N) =  binomial
// fL = ln(f(p))
// fL' = d(ln(f(p))/dp =  (k/p) - ( (n-k)/(1-p))
// fL'' = d2(ln(f(p))/d2p =  (-k/p^2)  - ( (n-k) / ((1-p)^2) )
// f' = df(p)/dp  =  f(p) * fL'
// f'' = d2f(p)/d2p = f(p) * (fL'^2 +fL'')

function rbinomDerivative($varray,$N,$m,$iProb=1,$iRate=0) {

  $derivs=[];
  $derivs2=[];

  $robj= rbinomPdf_init($N,$m);   // not needed if each row is [rate,prob]
  foreach ($varray as $ij=>$va1) {
      if (!is_array($va1) ) {                // scalar
          $arates=$va1;
          $aprob=rbinomPdf($robj,$va1);
      } else {
         $arate=$va1[$iRate];
         $aprob=$va1[$iProb];
     }
     if ($arate==0.0 || $arate==1.0) {  // edge cases
         $derivs[$ij]=[$arate,0,0];
         continue ;
    }
// now the FL' and FL''
     $t2a=($m/$arate); $t2b=($N-$m) / (1-$arate);
     $t2=$t2a-$t2b;
     $fLp=$t2;

     $t3a=(-$m/($arate**2)) ; $t3b=  ($N-$m)  / ( (1-$arate)**2) ;
     $fLpp=$t3a-$t3b;

     $aderiv=$aprob*$fLp ;
     $aderiv2=$aprob* (($fLp**2) + $fLpp)  ;
     $derivs[$ij]=[$arate,$aderiv,$aderiv2,$fLp,$fLpp];


   }
   return $derivs;
}
