# Modern Roman numerals

These two scripts, written in GNU Octave/Matlab language, convert Roman numerals to and from any integer (scalar, vector or matrix), including large numbers greater than 4999 with the parenthesis notation (multiplies by 1000).

## num2roman.m
The function `num2roman` uses strict rules of modern notation (substractive principle for 4 and 9 bases) except for the common `MMMM` form replacing `(IV)`.
### Examples
```matlab
>> num2roman(1968)
MCMLXVIII
>> % first power of 10
>> num2roman(10.^(0:7))
   'I'    'X'    'C'    'M'    '(X)'    '(C)'    '(M)'    '((X))'
>> % from 1 to 100
>> reshape(num2roman(1:100),10,10)'
'I'        'II'        'III'        'IV'        'V'        'VI'        'VII'        'VIII'        'IX'        'X'   
'XI'       'XII'       'XIII'       'XIV'       'XV'       'XVI'       'XVII'       'XVIII'       'XIX'       'XX'  
'XXI'      'XXII'      'XXIII'      'XXIV'      'XXV'      'XXVI'      'XXVII'      'XXVIII'      'XXIX'      'XXX'
'XXXI'     'XXXII'     'XXXIII'     'XXXIV'     'XXXV'     'XXXVI'     'XXXVII'     'XXXVIII'     'XXXIX'     'XL'  
'XLI'      'XLII'      'XLIII'      'XLIV'      'XLV'      'XLVI'      'XLVII'      'XLVIII'      'XLIX'      'L'   
'LI'       'LII'       'LIII'       'LIV'       'LV'       'LVI'       'LVII'       'LVIII'       'LIX'       'LX'  
'LXI'      'LXII'      'LXIII'      'LXIV'      'LXV'      'LXVI'      'LXVII'      'LXVIII'      'LXIX'      'LXX'
'LXXI'     'LXXII'     'LXXIII'     'LXXIV'     'LXXV'     'LXXVI'     'LXXVII'     'LXXVIII'     'LXXIX'     'LXXX'
'LXXXI'    'LXXXII'    'LXXXIII'    'LXXXIV'    'LXXXV'    'LXXXVI'    'LXXXVII'    'LXXXVIII'    'LXXXIX'    'XC'  
'XCI'      'XCII'      'XCIII'      'XCIV'      'XCV'      'XCVI'      'XCVII'      'XCVIII'      'XCIX'      'C'   

```
## roman2num.m
`roman2num` is more flexible and is able to convert some other Roman notation possibilities, for instance the 3 different expressions `IC`,`XCIX`, and `XCVIIII` all return `99`:
```matlab
>> roman2num({'IC','XCIX','XCVIIII'})
    99    99    99
```
or the 2 expressions `MDXV` and `MCCCCCXV` both return `1515`:
```matlab
>> roman2num({'MDXV','MCCCCCXV'})
   1515   1515
```

## Author
**François Beauducel**, [IPGP](www.ipgp.fr), [beaudu](https://github.com/beaudu), beauducel@ipgp.fr

## Documentation
Type `doc num2roman` and `doc roman2num` for help and syntax.
