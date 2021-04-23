<?php
namespace wSurvey\Math\sortSearch ;


// sort and seearch (array and matrix) functions (wSurvey php library).

//============
// quick merge of two matrices, using a column to determine merging.
// Both must be sorted in increasing order
// Returned matrix is combination of gvals and gadds, sort on column ith
// Note that gvals is often, but does NOT have to be, a longer matrix whose first value is lower than the first value of gadds.

function merge_ordered($gvals,$gadds,$ith=0) {
  $igAdd=0;
  $igExist=0;
  $igNew=0;
  $ggNew=[];

// trivial cases
   if (count($gvals)>0 && count($gadds)==0) return $gvals;
   if (count($gadds)>0 && count($gvals)==0) return $gadds;

   if ($gvals[0][$ith]<$gadds[0][$ith]) {
      $ggNew[0]=$gvals[0];      // start with first value in gvals
      $igExist=1;
   } else {            // gvals has smallest value
      $igNew=-1;
      for ($ig=0;$ig<count($gadds);$ig++) {
          if ($gadds[$ig][$ith]>$gvals[0][$ith]) break; // never true on ig=0
          $igNew++;
          $ggNew[$igNew]=$gadds[$ig];
          $igAdd=$ig;
      }
   }
  while (1) {                 // keep going until end of origijnal (gval) or addins (gadds). Then break
    if ($igExist>=count($gvals)) break ;

    $gg1=$gvals[$igExist][$ith];
    $gg2=$gadds[$igAdd][$ith];
    if ($gg1<=$gg2) {    // less than next add in
        $igNew++;
        $ggNew[$igNew]=$gvals[$igExist];
        $igExist++;
        continue;
    }               // else , add the next available addin
    $igNew++;
    $ggNew[$igNew]=$gadds[$igAdd];
    $igAdd++;
    if ($igAdd>= count($gadds) ) break ;   // all done, but check if anythign left in gvals
  }
  if ($igExist<count($gvals) ) {   // only one of these is true
     for ($ifoo=$igExist;$ifoo<count($gvals);$ifoo++) {
         $ggNew[]=$gvals[$ifoo];
     }
  }
  if ($igAdd<count($gadds))  {
     for ($ifoo=$igAdd;$ifoo<count($gadds);$ifoo++) {
         $ggNew[]=$gadds[$ifoo];
     }
  }

  return $ggNew;
}

//=====================
// binary search in a sorted matrix -- each row is a vector --
//  search on the jcol'th column (of each row. If jcol not specified, or <0) matrix is a simple array
// Limit search to Barray[Bleft] to Barrary[Bright] (all elements if these are not specified
// update of  binarySearchIndex
//  returns index that brackets a value (the closest value in array that is <= the target value)
// returns -1 if  (Bvalue<Barray[0], count($Barray)-1 if Bvalue>end(Barray)  ...
//       otherwise ii s.t. Barray[ii]<=Bvalue<Barray[ii+1]
// special casses:  bvalue=Barray[0] : return 0.  bvalue=Barray[last]  : return last index (=count(Barray)-1)
//  Barray : the matrix
//  Bvalue :  the value to find (the index of Barray st Barray[index]<=BValue<Barray[index+1]
//  $ith  : column of Barray to examine (defult =1)
//  $Bleft : start searching at Barray[$Bleft]. Default it 0 (first element)
//  $Bright  : stop searcing at Barray[$Bright]. Default,  is end (count($Barray)-1). If <0, this many from end

function binarySearchMatrix($Barray, $Bvalue,$jcol=-1,$Bleft=0,$Bright=-1 ) {

   if ($Bright<=0) $Bright=count($Barray)+$Bright;  // portion of array to search in
   $Bright=min($Bright,count($Barray)-1);
   if ($Bleft<0) return null ;
   if ($Bright<$Bleft) return null ; // an error

   if ($Bright<0) {return -1 ;  }  // empty array  - so  less than first value"
     $bleftVal= ($jcol<0) ? $Barray[$Bleft] :  $Barray[$Bleft][$jcol] ;
       if  ($bleftVal==$Bvalue) {return 0; }  //  equals first value
       if ($bleftVal>$Bvalue) {return -1; }  //  before first value
     $brightVal= ($jcol<0) ? $Barray[$Bright] :  $Barray[$Bright][$jcol] ;
       if ($brightVal<=$Bvalue) return $Bright  ;   // after, or egual, to last value

    if ($Bright==1) return 0;     // 2 element array, and ge first and lt last ..  must be in between (after first element at [0])

    while ($Bleft <= $Bright) {
      if ($Bleft+1>=$Bright) return $Bleft  ;  {       // indices are equal, or next to each other --so return bleft (lower end)
      }
      $Bmidpoint = (int) floor(($Bleft + $Bright) / 2);
      $vmid= ($jcol<0) ? $Barray[$Bmidpoint] : $Barray[$Bmidpoint][$jcol]  ;
      $jj=0;
      if ($vmid<$Bvalue) {
             $jj=-1;
      } else if  ($vmid>$Bvalue) {  ;
             $jj=1;
       }

      if ($jj==-1) {           // The midpoint value is less than the value.
        $Bleft = $Bmidpoint ;
      } elseif ($jj==1) {        // The midpoint value is greater than the value.
         $Bright = $Bmidpoint ;
      } else {                                      // exact match?
        return $Bmidpoint ;             // exact match at this element of barray (so this element is lower bound)
      }

     }
     return NULL;   // should never happen

}


//=====================
// binary search, with optional callback for custom comparison
// if $f is specified, it is called with 2 arguments:
//   $f(Barray[icheck],$Bvalue)      (check the value of the icheck index in Barray against the value of Bvalue)
// and must return: -1,0, or 1  --
//      -1 : Barray[icheck]<Bvalue, 0 if = , 1 if >
//  Note that it is up to the caller to be sure that the arguments types expected by $f are the same types provided as the arguments to binarySearch
// If $f is not specified < and > are used
// returns -1 if  (Bvalue<Barray[0], count($Barray)-1 if end(Barray)>Bvalue ... otherwise ii s.t. Barray[ii]<=Bvalue<Barray[ii+1]
function binarySearchIndex($Barray, $Bvalue,$f=null) {

    $Bleft = 0;
    $Bright = count($Barray) - 1;     // Set the right pointer to the length of the array -1.

    $gotFunc=0;
    if (!is_null($f) && is_callable($f))   $gotFunc=1;

    if ($Bright<0) {return -1 ;  }  // empty array  - so "less than first value"

    if ($gotFunc==1) {
        if ($f($Barray[0],$Bvalue)==0) return 0 ;   // = first value
        if ($f($Barray[0],$Bvalue)==1) return -1 ;   // before first value
        if ($f($Barray[$Bright],$Bvalue)<=0) return $Bright  ;   // after or equal to last value
    } else {

       if (($Barray[0])==$Bvalue) {return 0; }  //  before first value
       if (($Barray[0])>$Bvalue) {return -1; }  //  before first value
       if (($Barray[$Bright])<=$Bvalue) return $Bright  ;   // after, or egual, to last value
    }

    if ($Bright==1) return 0;     // 2 element array, and ge first and lt last ..  must be in between (after first element at [0])

    while ($Bleft <= $Bright) {
      if ($Bleft+1>=$Bright) {       // indices are equal, or next to each other
           return $Bleft  ; // in between bleft and bright, so return bleft (lower end)
      }
      $Bmidpoint = (int) floor(($Bleft + $Bright) / 2);
      if ($gotFunc!==1) {
          $vmid=$Barray[$Bmidpoint]  ;
          $jj=0;
          if ($vmid<$Bvalue) {
             $jj=-1;
          } else if  ($vmid>$Bvalue) {  ;
           $jj=1;
          }
      } else {          // usefunc
        $jj=$f($Barray[$Bmidpoint],$Bvalue);
      }
      if ($jj==-1) {           // The midpoint value is less than the value.
        $Bleft = $Bmidpoint ;
      } elseif ($jj==1) {        // The midpoint value is greater than the value.
         $Bright = $Bmidpoint ;
      } else {                                      // exact match?
        return $Bmidpoint ;
      }

     }
    return NULL;   // should never happen

}

//===========
// golden rule search adapted from https://en.wikipedia.org/wiki/Golden-section_search
//    Given a function f with a single local minimum in the interval [a,b]
//     grSearch returns a subset interval  that contains the minimum.
//     $f : function to minimize.  MUST be an anonymous function declared BEFORE this is called, that accepts one or two arguments
//          Or: a 2 element array: $f[0] = the anonymous function, $f[1] = other argument(s)
//     $a0 : low end of range   (x values)
//     $b0  : upper end of range
//     $xtol : stop if  the two "golden rule points"  are less that this distance apart. The subset interval will be a bit larger than this
//            If not specified, 1e-5 is used.
//            Special case: if $xtol<0, $xtol is the number of steps to take.
//
//  Returns:
//     [$x0,$x1,$y,$n,$stopHow]
//  x0 and x1 : the top and bottom of the range (within which the search stopped
//   y : the value at the stopping x
//   n : number of steps taken
//       stopHow:  stop condition code
//             0: a0 and b0 are close, so just analyze at a0 and return
//             1 : xtol met. Return value at x0
//             2 :  xtol met. Return value at x1
//             3 :  max (or specified) iterations. Return value at x0
//             4 :  max (or specified) iterations. Return value at x1

//   how: 1,2 3, or 4.     1 or 2 : tolerance met. 3 or 4: # iteratiosn met.  1 and : stop after computing lower end, 2 and 4 after upper.

function grSearch($f0,$a0,$b0,$xtol=1e-5) {
  $invphi = (sqrt(5) - 1) / 2   ; //   0.6180339
  $invphi2 = (3 - sqrt(5)) / 2  ; // 0.3819660
  $a=min($a0,$b0);   $b=max($a0,$b0);          // a is point minX (point 1)$b is point maxX (point 4)
  $h = $b - $a   ;
  if (is_array($f0)) {
      $f=$f0[0];
      $z=$f0[1];        // additional argument(s) to set to f0 function (these do not change)
      $do2=1;
  } else {
      $f=$f0;
      $do2=0;
  }
  if ($xtol<0) {          // special case: use set # of iterations
      $n=intVal(-$xtol);  $xtol=0;
  } else {
      $n = intVal(ceil(log($xtol / $h) / log($invphi))) ;   // 'Required' steps to achieve tolerance
  }
  if ($h <= $xtol) {     // odd case, bounds very close, so just analyze at lower bound and return
     if ($do2==0) {
         $y0=$f($a);
     } else {
       $y0=$f($a,$z)  ;
     }
     return [$a,$b,$y0,0,0] ;
  }

    $c = $a + $invphi2 * $h    ;    // point #2
    $d = $a + $invphi * $h    ;     // point #3

     if ($do2==0) {
        $yc = $f($c)     ;
        $yd = $f($d)     ;
     } else {
        $yc = $f($c,$z)     ;
        $yd = $f($d,$z)     ;
     }

    for ($k=0;$k<$n;$k++) {    // up to $n iterations
        if ($yc < $yd) {
            $b = $d  ;                      // max = old point3. minx does not change
            $d = $c ;                       // point 3 = old point 2
            $yd = $yc    ;
            $h = $invphi * $h   ;
            $c = $a + $invphi2 * $h    ;  //  point 2 is between minx and  new point3
            $yc = ($do2==0) ?   $f($c)  : $f($c,$z)  ;
        } else {
            $a = $c  ;    // minx = old point 2
            $c = $d  ;    // point2 = old point 3
            $yc = $yd   ;
            $h = $invphi * $h   ;
            $d = $a + $invphi * $h  ;   // point3 is between old point 2 and old (still used) max
            $yd = ($do2==0) ?   $f($d)  : $f($d,$z)  ;
        }

        if (($b-$a)<$xtol) {         // if using set $n, $xtol=0 so never satisfied
          if ($yc < $yd) {
              return [$a,$d,$yc,$k,1]  ;
          } else {
              return [$c,$b,$yd,$k,2] ;
          }
       }

    }
    if ($yc < $yd) {
        return [$a,$d,$yc,$n,3]  ;
    } else {
        return [$c,$b,$yd,$n,4] ;
    }

}


//==========
//https://php.net/manual/en/function.array-multisort.php
//  
//   $newArray=array_msort($array,$colSortInstructions,$sortType)
// 
//   where
//      $array: the associative array to sort
//      $colSortInstructions: an associative array congaining instructions on how to sort $array -- which columns to sort on first, in what order
//                              Each element in this array is  'colName'=>sortOrder
//                              sortOrder is either SORT_ASC or SORT_DESC (php constants, so do NOT preced with a $)
//      $sortType: optional. The type of sort (string, numeric, etc)
//                           This can be one of: SORT_REGULAR, SORT_NUMERICl, SORT_STRING, SORT_LOCALE_STRING, SORT_NATURAL, SORT_FLAG_CASE (acan be OR with SORT_STRING or SORT_NATURAL)
//                           If not specified, SORT_REGULAR is used
//
//  Returns sorted array. Thus, $array is NOT changed.
//
//  Example (assumes $arr1 had an 'id1' and 'id2' colum, perhaps other columns also
//       $newarrray=array_msort($arr1, array('id1'=>SORT_DESC, 'id2'=>SORT_ASC'],SORT_STRING);
//
//  Sorted in order of indices. So first sort on "id1", and then within blocks of data with the same value of 'id1', sort by 'id2'
//
// See php array_multisort for descriptiong of the sortOrder and sortType options

function array_msort($array, $cols,$sortType=SORT_REGULAR,$debug=0) {
    $colarr = array();
    foreach ($cols as $col => $order) {
        $colarr[$col] = array();
//        foreach ($array as $k => $row) { $colarr[$col]['_'.$k] = strtolower($row[$col]); }
        foreach ($array as $k => $row) {
            $colarr[$col]['_'.$k] = $row[$col] ; 
        }
    }
    $eval = 'array_multisort(';
    foreach ($cols as $col => $order) {
        $eval .= '$colarr[\''.$col.'\'],'.$order.','.$sortType.',';
    }
    $eval = substr($eval,0,-1).');';

    if ($debug== 1) {
      print '<BR>SORT_ASC='.SORT_ASC.', SORT_DESC= '.SORT_DESC.' ,  SORT_NUMERIC ='.SORT_NUMERIC.', SORT_REGULAR= '.SORT_REGULAR.', SORT_STRING='.SORT_STRING ;
      print '<div style="overflow:scroll;xheight:14em">';
      print '<br>original ';
      do_dump($array);
      print '<br> sort instructions ';
      do_dump($colarr);
      print '<br>the call ';
      do_dump($eval);
     }

    eval($eval);
    $ret = array();
    foreach ($colarr as $col => $arr) {
        foreach ($arr as $k => $v) {
            $k = substr($k,1);
            if (!isset($ret[$k])) $ret[$k] = $array[$k];
            $ret[$k][$col] = $array[$k][$col];
        }
    }

    if ($debug== 1) {
       print '<br>Sorted ';
       do_dump($ret);
       print '</div>';
       exit;
    }

    return $ret;

}

