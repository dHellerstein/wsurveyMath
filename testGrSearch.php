<?php
// Demo of grSearch (golden rule search)

require('sortSearch.php' );    // or whereer you chose to store the  wSurvey\Math release.

use wSurvey\Math\sortSearch  as ss ;

$coeffs=['B0'=>100,'B1'=>-25,'B2'=>3];  // 2nd degree polynomial coefficients.  Note that min will be at -25 + 6*x = 0  == 4.166666

// the anonymous function
$fPoly=function($aval,$acoeffs) {
  $v1=$acoeffs['B0'] + ($acoeffs['B1']*$aval) + ($acoeffs['B2']*($aval**2)) ;
  return $v1;
}    ;

// -- version used to display a few values
$fPoly2=function($aval) {
$bcoeffs=['B0'=>100,'B1'=>-25,'B2'=>3];  // 2nd degree polynomial coefficients
  $v1=$bcoeffs['B0'] + ($bcoeffs['B1']*$aval) + ($bcoeffs['B2']*($aval**2)) ;
  return $v1;
}    ;

$tst=[0.,.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5];

$goo=array_map($fPoly2,$tst);

 echo "<br> golden rule results of polynomial: ".$coeffs['B0'] . ' + ' . $coeffs['B1'].' * x +  '. $coeffs['B2']. ' * x^2 ';

echo "<br>Some values of this : ";
echo '<ul>';
foreach ($tst as $j=>$av) echo '<li> f('.$av.') = ' .$goo[$j];
echo '</ul>';

$res=ss\grSearch([$fPoly,$coeffs],0,10000)   ;     //  use default tolerance
echo "<br>Search results [xlower, xupper, y, #steps, stopCondition]: ".implode(', ',$res);


?>
