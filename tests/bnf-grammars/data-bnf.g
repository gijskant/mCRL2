DataExpr
  : [da_id] Id
  | [da_numb] Number
  | [da_true] 'true'
  | [da_false] 'false'
 // | '[' ']'
 // | '{' '}'
 // | '[' DataExprList ']'
 // | '{' BagEnumEltList '}'
//  | '{' VarDecl '|' DataExpr '}'
//  | '{' DataExprList '}'
  | [da_brk] '(' DataExpr ')'
//  | DataExpr '[' DataExpr '->' DataExpr ']'  $unary_left  13
//  | DataExpr '(' DataExprList ')'            $unary_left  13
  | [da_neg] '!' DataExpr                             $unary_right 12
  | [da_inv] '-' DataExpr                             $unary_right 12
 // | '#' DataExpr                             $unary_right 12
  | [da_forall] 'forall' VarsDeclList '.' DataExpr    
  | [da_exists] 'exists' VarsDeclList '.' DataExpr     
//  | [da_lamb] 'lambda' VarsDeclList '.' DataExpr    
  | [da_impl] DataExpr '=>' DataExpr
  | [da_disj] DataExpr '||' DataExpr
  | [da_conj] DataExpr '&&' DataExpr
  | [da_eq] DataExpr '==' DataExpr
  | [da_neq] DataExpr '!=' DataExpr
  | [da_lt] DataExpr '<' DataExpr
  | [da_leq] DataExpr '<=' DataExpr
  | [da_geq] DataExpr '>=' DataExpr
  | [da_gt] DataExpr '>' DataExpr
 // | DataExpr ('in' $binary_op_left 6) DataExpr
 // | DataExpr ('|>' $binary_op_right 7) DataExpr
//  | DataExpr ('<|' $binary_op_left 8) DataExpr
//  | DataExpr ('++' $binary_op_left 9) DataExpr
  | [da_plus] DataExpr '+' DataExpr
  | [da_min] DataExpr '-' DataExpr
  | [da_divides] DataExpr '/' DataExpr
 // | DataExpr ('div' $binary_op_left 11) DataExpr
 // | DataExpr ('mod' $binary_op_left 11) DataExpr
  | [da_times] DataExpr '*' DataExpr
 // | DataExpr ('.' $binary_op_left 12) DataExpr
//  | [da_whr] DataExpr 'whr' AssignmentList 'end'      $unary_left 0
  ;
  
// Needed to get arguments for prop var that is 
DataExprListInst: [dalinst] DataExpr
                | [dallinst] DataExpr ',' DataExprListInst
			    ;