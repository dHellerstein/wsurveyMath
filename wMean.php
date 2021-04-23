<?php
namespace wSurvey\Math\wMean ;

// weighted mean and variance, with selection
//   returns [mean, variance, #of obs used, cumulative weight]
//  x1: array (any style of indices) to compute statistics (mean, sd, #obs) on
//      Non numeric values are dropped.
//  weights: weights used when computing stats.
//     If not specified, no weighting.
//      If [], no weighting.
//      If  index (into x1) does not exist in weights, or is non numeric,use a weight of 1
//
//  Note that weights are "normalize" -- so their sum adds to 1 (i.e.; the weights are a discrete probablity distribution function)
//  selx : selection values (array)
//  selVal : selection criteria (scaler)
//     If selx or selval not specified, all values in x1 are used
//     If selx=[], all valeus in x1 are used
//     If an index (into x1) does not exists in selx, the value is used
//     Otherwise, if selx[i]<>selval, the observation is dropped (not used)
//     Not used observations are NOT included in the "# of obs used" returned value
//
// Tecnhical notes:
//  *  when computing means and vairiabes weighting means "replicating rows" --
//     so a weight of 2 implies "there are two observations identical to this one"
//     Fractional weights are allowed (such  as 3.561)
//  * Since this creates internal copies of x1, it is not particularly efficient.


 function wcmp_mean_var($x1,$weights=[],$selx=[],$selVal='' ) {
   $nx0=count($x1);
   if ($nx0==0) return [0,0,0,0] ; //0 rows means unable to calculate
   $xuse=[]; $wuse=[];

   $nselx=count($selx);
   if (func_num_args()<4) $nselx=0;

   $nweight=count($weights);
   
// create internal xuse and wuse arrays (dropping obs, filling in weights if not specified)
   foreach ($x1 as $ix1=>$ax) {
      if (!is_numeric($ax)) continue ;
      $wuse1=1.0;
      if ($nselx>0) {            // check selx?
         if (!array_key_exists($ix1,$selx)) {   // keep if no match in selx
            if ($nweight>0) {
               if (array_key_exists($ix1,$weights)) $wuse1=floatval($weights[$ix1]);
               $wuse[]=$wuse1;
            }
            $wuse[]=$wuse1;
            $xuse[]=$ax;
            continue;
         }
         if ($selx[$ix1]!=$selVal) continue;   // do not retain
         if ($nweight>0) {
            if (array_key_exists($ix1,$weights)) $wuse1=floatval($weights[$ix1]);
            $wuse[]=$wuse1;
         }
         $xuse[]=floatVal($ax);;
     }   else {     // nselx =0 -- keep!
         if ($nweight>0) {
            if (array_key_exists($ix1,$weights)) $wuse1=floatval($weights[$ix1]);
            $wuse[]=$wuse1;
         }
         $xuse[]=floatVal($ax);
     }
   }    // for


   $nx=count($xuse);          // after removal of obs
   if ($nx==0) return [0,0,0,0] ; //0 rows means unable to calculate
   if ($nx==1) return [$xuse[0],0,1,1] ; //0 rows means unable to calculate  variance

// if no weighting, can use shortcuts     http://www.statisticslectures.com/topics/variancesample/
   if ($nweight==0) {
     $x1Sum=0; $x1SumSquare=0;
     for ($ix=0;$ix<$nx;$ix++) {
       $ax1=$xuse[$ix];
       $x1Sum+=$ax1 ; $x1SumSquare+=$ax1*$ax1 ;
     }
     $num1b=($x1Sum*$x1Sum);
     $xmean=$x1Sum/$nx;
     $var1=($x1SumSquare-($num1b/$nx))/($nx) ;  // could use nx-1 (sample mean), but for compatabilioty with unweighted, don't
     return [$xmean,$var1,count($xuse),0];
   }

// weighting!

   $x1Sum=0; $x1SumSquare=0; $totCt=0;
   for ($ix=0;$ix<$nx;$ix++) {
       $ax1=floatval($xuse[$ix]);
       $aw1=floatval($wuse[$ix]);
       $totCt+=$aw1;             // this is used to "normalize" the weights
       $x1Sum+=($ax1*$aw1) ;
    }
    $wmean=$x1Sum/$totCt ;    // note use of totCt to normalize the weights

    $sumXdiffSq=0;
    $aw1=1.0;
   for ($ix=0;$ix<$nx;$ix++) {
       $ax1=$xuse[$ix];
       if ($nweight>0) $aw1=$wuse[$ix];
       $xdiffSq=pow($ax1-$wmean,2);
       $sumXdiffSq+=($xdiffSq*$aw1);
   }
   $wvar=$sumXdiffSq/$totCt;
   return [$wmean,$wvar,$nx,$totCt];
}


