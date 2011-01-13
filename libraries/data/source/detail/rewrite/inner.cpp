// Author(s): Muck van Weerdenburg
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file inner.cpp

#include "boost.hpp" // precompiled headers

#include "mcrl2/data/detail/rewrite/inner.h"

#define NAME "rewr_inner"

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <stdexcept>
#include <memory.h>
#include "mcrl2/atermpp/aterm_access.h"
#include "mcrl2/core/detail/memory_utility.h"
#include "mcrl2/core/detail/struct_core.h"
#include "mcrl2/core/print.h"
#include "mcrl2/core/messaging.h"
#include "mcrl2/core/aterm_ext.h"
#include "mcrl2/data/data_equation.h"

using namespace mcrl2::core;
using namespace mcrl2::core::detail;

namespace mcrl2 {
  namespace data {
    namespace detail {

static AFun nilAFun;
static AFun opidAFun;
static AFun ruleAFun;
static unsigned int is_initialised = 0;

static AFun afunS, afunM, afunF, afunN, afunD, afunR, afunCR, afunC, afunX, afunRe, afunCRe, afunMe;
static ATerm dummy;

#define isS(x) ATisEqualAFun(ATgetAFun(x),afunS)
#define isM(x) ATisEqualAFun(ATgetAFun(x),afunM)
#define isF(x) ATisEqualAFun(ATgetAFun(x),afunF)
#define isN(x) ATisEqualAFun(ATgetAFun(x),afunN)
#define isD(x) ATisEqualAFun(ATgetAFun(x),afunD)
#define isR(x) ATisEqualAFun(ATgetAFun(x),afunR)
#define isCR(x) ATisEqualAFun(ATgetAFun(x),afunCR)
#define isC(x) ATisEqualAFun(ATgetAFun(x),afunC)
#define isX(x) ATisEqualAFun(ATgetAFun(x),afunX)
#define isRe(x) ATisEqualAFun(ATgetAFun(x),afunRe)
#define isCRe(x) ATisEqualAFun(ATgetAFun(x),afunCRe)
#define isMe(x) ATisEqualAFun(ATgetAFun(x),afunMe)

#define ATAgetArgument(x,y) ((ATermAppl) ATgetArgument(x,y))
#define ATLgetArgument(x,y) ((ATermList) ATgetArgument(x,y))
#define ATAgetFirst(x) ((ATermAppl) ATgetFirst(x))
#define ATLgetFirst(x) ((ATermList) ATgetFirst(x))
#define gsIsOpId(x) (ATgetAFun(x) == opidAFun)

#define ATisList(x) (ATgetType(x) == AT_LIST)
#define ATisAppl(x) (ATgetType(x) == AT_APPL)
#define ATisInt(x) (ATgetType(x) == AT_INT)

#define is_nil(x) (ATisList(x)?false:(ATgetAFun((ATermAppl) x) == nilAFun))


static void initialise_common()
{
  if ( is_initialised == 0 )
  {
    nilAFun = ATgetAFun(gsMakeNil());
    ATprotectAFun(nilAFun);
    opidAFun = ATgetAFun(static_cast<ATermAppl>(sort_bool::true_()));
    ATprotectAFun(opidAFun);
    ruleAFun = ATmakeAFun("@RULE@",4,ATfalse);
    ATprotectAFun(ruleAFun);

    afunS = ATmakeAFun("@@S",2,ATfalse); // Store term ( target_variable, result_tree )
    ATprotectAFun(afunS);
    afunM = ATmakeAFun("@@M",3,ATfalse); // Match term ( match_variable, true_tree , false_tree )
    ATprotectAFun(afunM);
    afunF = ATmakeAFun("@@F",3,ATfalse); // Match function ( match_function, true_tree, false_tree )
    ATprotectAFun(afunF);
    afunN = ATmakeAFun("@@N",1,ATfalse); // Go to next parameter ( result_tree )
    ATprotectAFun(afunN);
    afunD = ATmakeAFun("@@D",1,ATfalse); // Go down a level ( result_tree )
    ATprotectAFun(afunD);
    afunR = ATmakeAFun("@@R",1,ATfalse); // End of tree ( matching_rule )
    ATprotectAFun(afunR);
    afunCR = ATmakeAFun("@@CR",2,ATfalse); // End of tree ( condition, matching_rule )
    ATprotectAFun(afunCR);
    afunC = ATmakeAFun("@@C",3,ATfalse); // Check condition ( condition, true_tree, false_tree )
    ATprotectAFun(afunC);
    afunX = ATmakeAFun("@@X",0,ATfalse); // End of tree
    ATprotectAFun(afunX);
    afunRe = ATmakeAFun("@@Re",2,ATfalse); // End of tree ( matching_rule , vars_of_rule)
    ATprotectAFun(afunRe);
    afunCRe = ATmakeAFun("@@CRe",4,ATfalse); // End of tree ( condition, matching_rule, vars_of_condition, vars_of_rule )
    ATprotectAFun(afunCRe);
    afunMe = ATmakeAFun("@@Me",2,ATfalse); // Match term ( match_variable, variable_index )
    ATprotectAFun(afunMe);
    dummy = NULL;
    ATprotect(&dummy);
    dummy = (ATerm) gsMakeNil();
  }

  is_initialised++;
}

static void finalise_common()
{
  assert(is_initialised > 0 );
  is_initialised--;

  if ( is_initialised == 0 )
  {
    ATunprotect(&dummy);
    ATunprotectAFun(afunMe);
    ATunprotectAFun(afunCRe);
    ATunprotectAFun(afunRe);
    ATunprotectAFun(afunX);
    ATunprotectAFun(afunC);
    ATunprotectAFun(afunCR);
    ATunprotectAFun(afunR);
    ATunprotectAFun(afunD);
    ATunprotectAFun(afunN);
    ATunprotectAFun(afunF);
    ATunprotectAFun(afunM);
    ATunprotectAFun(afunS);

    ATunprotectAFun(ruleAFun);
    ATunprotectAFun(opidAFun);
    ATunprotectAFun(nilAFun);
  }
}


#ifdef _INNER_STORE_TREES
size_t RewriterInnermost::write_tree(FILE *f, ATermAppl tree, size_t *num_states)
{
  if ( isS(tree) )
  {
    size_t n = write_tree(f,ATAgetArgument(tree,1),num_states);
    fprintf(f,"n%i [label=\"S(%s)\"]\n",*num_states,ATgetName(ATgetAFun(ATAgetArgument(ATAgetArgument(tree,0),0))));
    fprintf(f,"n%i -> n%i\n",*num_states,n);
    return (*num_states)++;
  } else if ( isM(tree) )
  {
    size_t n = write_tree(f,ATAgetArgument(tree,1),num_states);
    size_t m = write_tree(f,ATAgetArgument(tree,2),num_states);
    if ( ATisInt(ATgetArgument(tree,0)) )
    {
      fprintf(f,"n%i [label=\"M(%i)\"]\n",*num_states,ATgetInt((ATermInt) ATgetArgument(tree,0)));
    } else {
      fprintf(f,"n%i [label=\"M(%s)\"]\n",*num_states,ATgetName(ATgetAFun(ATAgetArgument(ATAgetArgument(tree,0),0))));
    }
    fprintf(f,"n%i -> n%i [label=\"true\"]\n",*num_states,n);
    fprintf(f,"n%i -> n%i [label=\"false\"]\n",*num_states,m);
    return (*num_states)++;
  } else if ( isF(tree) )
  {
    size_t n = write_tree(f,ATAgetArgument(tree,1),num_states);
    size_t m = write_tree(f,ATAgetArgument(tree,2),num_states);
    if ( ATisInt(ATgetArgument(tree,0)) )
    {
      fprintf(f,"n%i [label=\"F(%s)\"]\n",*num_states,ATgetName(ATgetAFun(ATAgetArgument(int2term[ATgetInt((ATermInt) ATgetArgument(tree,0))],0))));
    } else {
      fprintf(f,"n%i [label=\"F(%s)\"]\n",*num_states,ATgetName(ATgetAFun(ATAgetArgument(ATAgetArgument(tree,0),0))));
    }
    fprintf(f,"n%i -> n%i [label=\"true\"]\n",*num_states,n);
    fprintf(f,"n%i -> n%i [label=\"false\"]\n",*num_states,m);
    return (*num_states)++;
  } else if ( isD(tree) )
  {
    size_t n = write_tree(f,ATAgetArgument(tree,0),num_states);
    fprintf(f,"n%i [label=\"D\"]\n",*num_states);
    fprintf(f,"n%i -> n%i\n",*num_states,n);
    return (*num_states)++;
  } else if ( isN(tree) )
  {
    size_t n = write_tree(f,ATAgetArgument(tree,0),num_states);
    fprintf(f,"n%i [label=\"N\"]\n",*num_states);
    fprintf(f,"n%i -> n%i\n",*num_states,n);
    return (*num_states)++;
  } else if ( isC(tree) )
  {
    size_t n = write_tree(f,ATAgetArgument(tree,1),num_states);
    size_t m = write_tree(f,ATAgetArgument(tree,2),num_states);
    gsfprintf(f,"n%i [label=\"C(%P)\"]\n",*num_states,fromInner(ATgetArgument(tree,0)));
    fprintf(f,"n%i -> n%i [label=\"true\"]\n",*num_states,n);
    fprintf(f,"n%i -> n%i [label=\"false\"]\n",*num_states,m);
    return (*num_states)++;
  } else if ( isR(tree) )
  {
    gsfprintf(f,"n%i [label=\"R(%P)\"]\n",*num_states,fromInner(ATgetArgument(tree,0)));
    return (*num_states)++;
  } else if ( isX(tree) )
  {
    ATfprintf(f,"n%i [label=\"X\"]\n",*num_states);
    return (*num_states)++;
  }

  return -1;
}

void RewriterInnermost::tree2dot(ATermAppl tree, char *name, char *filename)
{
  FILE *f;
  size_t num_states = 0;

  if ( (f = fopen(filename,"w")) == NULL )
  {
    perror("fopen");
    return;
  }

  fprintf(f,"digraph \"%s\" {\n",name);
  write_tree(f,tree,&num_states);
  fprintf(f,"}\n");

  fclose(f);
}
#endif

static void term2seq(ATerm t, ATermList *s, size_t *var_cnt)
{
  if ( ATisList(t) )
  {
    ATermList l;

    l = ATgetNext((ATermList) t);
    t = ATgetFirst((ATermList) t);

    *s = ATinsert(*s, (ATerm) ATmakeAppl3(afunF,(ATerm) t,dummy,dummy));

    for (; !ATisEmpty(l); l=ATgetNext(l))
    {
      term2seq(ATgetFirst(l),s,var_cnt);
      if ( !ATisEmpty(ATgetNext(l)) )
      {
        *s = ATinsert(*s, (ATerm) ATmakeAppl1(afunN,dummy));
      }
    }
    *s = ATinsert(*s, (ATerm) ATmakeAppl1(afunD,dummy));
  } else if ( ATisInt(t) )
  {
    term2seq((ATerm) ATmakeList1(t),s,var_cnt);
  } else if ( gsIsDataVarId((ATermAppl) t) )
  {
    ATerm store = (ATerm) ATmakeAppl2(afunS,(ATerm) t,dummy);

    if ( ATindexOf(*s,store,0) != ATERM_NON_EXISTING_POSITION)
    {
      *s = ATinsert(*s, (ATerm) ATmakeAppl3(afunM,(ATerm) t,dummy,dummy));
    } else {
      (*var_cnt)++;
      *s = ATinsert(*s, store);
    }
  } else {
    assert(0);
  }

}

static void get_used_vars_aux(ATerm t, ATermList *vars)
{
  if ( ATisList(t) )
  {
    for (; !ATisEmpty((ATermList) t); t=(ATerm) ATgetNext((ATermList) t))
    {
      get_used_vars_aux(ATgetFirst((ATermList) t),vars);
    }
  } else if ( ATisAppl(t) )
  {
    if ( gsIsDataVarId((ATermAppl) t) )
    {
      if ( ATindexOf(*vars,t,0) == ATERM_NON_EXISTING_POSITION )
      {
        *vars = ATinsert(*vars,t);
      }
    } else {
      size_t a = ATgetArity(ATgetAFun((ATermAppl) t));
      for (size_t i=0; i<a; i++)
      {
        get_used_vars_aux(ATgetArgument((ATermAppl) t,i),vars);
      }
    }
  }
}

static ATermList get_used_vars(ATerm t)
{
  ATermList l = ATmakeList0();

  get_used_vars_aux(t,&l);

  return l;
}

static ATermList create_sequence(ATermAppl rule, size_t *var_cnt, ATermInt trueint)
{
  ATermAppl pat = (ATermAppl) ATgetArgument(rule,2);
  ATerm cond = ATgetArgument(rule,1);
  ATerm rslt = ATgetArgument(rule,3);
  ATermList pars = ATmakeList0();
  ATermList rseq = ATmakeList0();

  pars = (ATermList) pat;
  //ATfprintf(stderr,"pattern pars: %t\n",pars);
  for (; !ATisEmpty(pars); pars=ATgetNext(pars))
  {
    term2seq(ATgetFirst(pars),&rseq,var_cnt);
    if ( !ATisEmpty(ATgetNext(pars)) )
    {
      rseq = ATinsert(rseq, (ATerm) ATmakeAppl1(afunN,dummy));
    }
  }
  //ATfprintf(stderr,"rseq: %t\n",rseq);
  if ( ATisInt(cond)/* && gsIsNil((ATermAppl) cond) */ && ATisEqual(cond, trueint)) // JK 15/10/2009 default condition is true
    rseq = ATinsert(rseq,(ATerm) ATmakeAppl2(afunRe,rslt,(ATerm) get_used_vars(rslt)));
  else
    rseq = ATinsert(rseq,(ATerm) ATmakeAppl4(afunCRe,cond,rslt,(ATerm) get_used_vars(cond),(ATerm) get_used_vars(rslt)));

  return ATreverse(rseq);
}


// Structure for build_tree paramters
typedef struct {
  ATermList Flist;   // List of sequences of which the first action is an F
  ATermList Slist;   // List of sequences of which the first action is an S
  ATermList Mlist;   // List of sequences of which the first action is an M
  ATermList stack;   // Stack to maintain the sequences that do not have to
                     // do anything in the current term
  ATermList upstack; // List of sequences that have done an F at the current
                     // level
} build_pars;

static void initialise_build_pars(build_pars *p)
{
  p->Flist = ATmakeList0();
  p->Slist = ATmakeList0();
  p->Mlist = ATmakeList0();
  p->stack = ATmakeList1((ATerm) ATmakeList0());
  p->upstack = ATmakeList0();
  /*ATprotectList(&p->Flist);
  ATprotectList(&p->Slist);
  ATprotectList(&p->Mlist);
  ATprotectList(&p->stack);
  ATprotectList(&p->upstack);*/
}

static void finalise_build_pars(build_pars * /*p*/)
{
/*	ATunprotectList(&p->Flist);
  ATunprotectList(&p->Slist);
  ATunprotectList(&p->Mlist);
  ATunprotectList(&p->stack);
  ATunprotectList(&p->upstack);*/
}

static ATermList add_to_stack(ATermList stack, ATermList seqs, ATermAppl *r, ATermList *cr)
{
  if ( ATisEmpty(stack) )
  {
    return stack;
  }

  ATermList l = ATmakeList0();
  ATermList h = ATLgetFirst(stack);

  for (; !ATisEmpty(seqs); seqs=ATgetNext(seqs))
  {
    ATermList e = ATLgetFirst(seqs);

    if ( isD(ATAgetFirst(e)) )
    {
      l = ATinsert(l,(ATerm) ATgetNext(e));
    } else if ( isN(ATAgetFirst(e)) )
    {
      h = ATinsert(h,(ATerm) ATgetNext(e));
    } else if ( isRe(ATAgetFirst(e)) )
    {
      *r = ATAgetFirst(e);
    } else {
      *cr = ATinsert(*cr,ATgetFirst(e));
    }
  }

  return ATinsert(add_to_stack(ATgetNext(stack),l,r,cr),(ATerm) h);
}

static void add_to_build_pars(build_pars *pars,ATermList seqs, ATermAppl *r, ATermList *cr)
{
  ATermList l = ATmakeList0();

  for (; !ATisEmpty(seqs); seqs=ATgetNext(seqs))
  {
    ATermList e = ATLgetFirst(seqs);

    if ( isD(ATAgetFirst(e)) || isN(ATAgetFirst(e)) )
    {
      l = ATinsert(l,(ATerm) e);
    } else if ( isS(ATAgetFirst(e)) )
    {
      pars->Slist = ATinsert(pars->Slist,(ATerm) e);
    } else if ( isMe(ATAgetFirst(e)) ) // M should not appear at the head of a seq
    {
      pars->Mlist = ATinsert(pars->Mlist,(ATerm) e);
    } else if ( isF(ATAgetFirst(e)) )
    {
      pars->Flist = ATinsert(pars->Flist,(ATerm) e);
    } else if ( isRe(ATAgetFirst(e)) )
    {
      *r = ATAgetFirst(e);
    } else {
      *cr = ATinsert(*cr,ATgetFirst(e));
    }
  }

  pars->stack = add_to_stack(pars->stack,l,r,cr);
}

static char tree_var_str[20];
static ATermAppl createFreshVar(ATermAppl sort,size_t *i)
{
  sprintf(tree_var_str,"var-%lu",(*i)++);
  return gsMakeDataVarId(gsString2ATermAppl(tree_var_str),sort);
}

static ATermList subst_var(ATermList l, ATermAppl old, ATerm new_term, ATerm num, ATermList substs)
{
  if ( ATisEmpty(l) )
  {
    return l;
  }

  ATermAppl head = (ATermAppl) ATgetFirst(l);
  l = ATgetNext(l);

  if ( isM(head) )
  {
    if ( ATisEqual(ATgetArgument(head,0),old) )
    {
      head = ATmakeAppl2(afunMe,new_term,num);
    }
  } else if ( isCRe(head) )
  {
    ATermList l = (ATermList) ATgetArgument(head,2);
    ATermList m = ATmakeList0();
    for (; !ATisEmpty(l); l=ATgetNext(l))
    {
      if ( ATisEqual(ATgetFirst(l),old) )
      {
        m = ATinsert(m,num);
      } else {
        m = ATinsert(m,ATgetFirst(l));
      }
    }
    l = (ATermList) ATgetArgument(head,3);
    ATermList n = ATmakeList0();
    for (; !ATisEmpty(l); l=ATgetNext(l))
    {
      if ( ATisEqual(ATgetFirst(l),old) )
      {
        n = ATinsert(n,num);
      } else {
        n = ATinsert(n,ATgetFirst(l));
      }
    }
    head = ATmakeAppl4(afunCRe,gsSubstValues(substs,ATgetArgument(head,0),true),gsSubstValues(substs,ATgetArgument(head,1),true),(ATerm) m, (ATerm) n);
  } else if ( isRe(head) )
  {
    ATermList l = (ATermList) ATgetArgument(head,1);
    ATermList m = ATmakeList0();
    for (; !ATisEmpty(l); l=ATgetNext(l))
    {
      if ( ATisEqual(ATgetFirst(l),old) )
      {
        m = ATinsert(m,num);
      } else {
        m = ATinsert(m,ATgetFirst(l));
      }
    }
    head = ATmakeAppl2(afunRe,gsSubstValues(substs,ATgetArgument(head,0),true),(ATerm) m);
  }

  return ATinsert(subst_var(l,old,new_term,num,substs),(ATerm) head);
}

//#define BT_DEBUG
#ifdef BT_DEBUG
#define print_return(x,y) ATermAppl a = y; ATfprintf(stderr,x "return %t\n\n",a); return a;
#else
#define print_return(x,y) return y;
#endif
//static size_t max_tree_vars;
static size_t *treevars_usedcnt;

static void inc_usedcnt(ATermList l)
{
  for (; !ATisEmpty(l); l=ATgetNext(l))
  {
    treevars_usedcnt[ATgetInt((ATermInt) ATgetFirst(l))]++;
  }
}

static ATermAppl build_tree(build_pars pars, size_t i)
{
#ifdef BT_DEBUG
ATfprintf(stderr,"build_tree(  %t  ,  %t  ,  %t  ,  %t  ,  %t  ,  %i  )\n\n",pars.Flist,pars.Slist,pars.Mlist,pars.stack,pars.upstack,i);
#endif

  if ( !ATisEmpty(pars.Slist) )
  {
    ATermList l,m;

    size_t k = i;
    ATermAppl v = createFreshVar(ATAgetArgument(ATAgetArgument(ATAgetFirst(ATLgetFirst(pars.Slist)),0),1),&i);
    treevars_usedcnt[k] = 0;

    l = ATmakeList0();
    m = ATmakeList0();
    for (; !ATisEmpty(pars.Slist); pars.Slist=ATgetNext(pars.Slist))
    {
      ATermList e = ATLgetFirst(pars.Slist);
       assert(k<((size_t)1)<<(8*sizeof(int)-1));
      e = subst_var(e,ATAgetArgument(ATAgetFirst(e),0),(ATerm) v,(ATerm) ATmakeInt((int)k),ATmakeList1((ATerm) gsMakeSubst(ATgetArgument(ATAgetFirst(e),0),(ATerm) v)));

      l = ATinsert(l,ATgetFirst(e));
      m = ATinsert(m,(ATerm) ATgetNext(e));
    }

    ATermAppl r = NULL;
    ATermList readies = ATmakeList0();

    pars.stack = add_to_stack(pars.stack,m,&r,&readies);

    if ( r == NULL )
    {
      ATermAppl tree;

      tree = build_tree(pars,i);
      for (; !ATisEmpty(readies); readies=ATgetNext(readies))
      {
        inc_usedcnt((ATermList) ATgetArgument(ATAgetFirst(readies),2));
        inc_usedcnt((ATermList) ATgetArgument(ATAgetFirst(readies),3));
        tree = ATmakeAppl3(afunC,ATgetArgument(ATAgetFirst(readies),0),(ATerm) ATmakeAppl1(afunR,ATgetArgument(ATAgetFirst(readies),1)),(ATerm) tree);
      }
      r = tree;
    } else {
      inc_usedcnt((ATermList) ATgetArgument(r,1));
      r = ATmakeAppl1(afunR,ATgetArgument(r,0));
    }

    if ( (treevars_usedcnt[k] > 0) || ((k == 0) && isR(r)) )
    {
      print_return("",ATmakeAppl2(afunS,(ATerm) v,(ATerm) r));
    } else {
      print_return("",r);
    }
  } else if ( !ATisEmpty(pars.Mlist) )
  {
    ATerm M = ATgetFirst(ATLgetFirst(pars.Mlist));

    ATermList l = ATmakeList0();
    ATermList m = ATmakeList0();
    for (; !ATisEmpty(pars.Mlist); pars.Mlist=ATgetNext(pars.Mlist))
    {
      if ( ATisEqual(M,ATgetFirst(ATLgetFirst(pars.Mlist))) )
      {
        l = ATinsert(l,(ATerm) ATgetNext(ATLgetFirst(pars.Mlist)));
      } else {
        m = ATinsert(m,ATgetFirst(pars.Mlist));
      }
    }
    pars.Mlist = m;

    ATermAppl true_tree,false_tree;
    ATermAppl r = NULL;
    ATermList readies = ATmakeList0();

    ATermList newstack = add_to_stack(pars.stack,l,&r,&readies);

    false_tree = build_tree(pars,i);

    if  ( r == NULL )
    {
      pars.stack = newstack;
      true_tree = build_tree(pars,i);
      for (; !ATisEmpty(readies); readies=ATgetNext(readies))
      {
        inc_usedcnt((ATermList) ATgetArgument(ATAgetFirst(readies),2));
        inc_usedcnt((ATermList) ATgetArgument(ATAgetFirst(readies),3));
        true_tree = ATmakeAppl3(afunC,ATgetArgument(ATAgetFirst(readies),0),(ATerm) ATmakeAppl1(afunR,ATgetArgument(ATAgetFirst(readies),1)),(ATerm) true_tree);
      }
    } else {
      inc_usedcnt((ATermList) ATgetArgument(r,1));
      true_tree = ATmakeAppl1(afunR,ATgetArgument(r,0));
    }

    if ( ATisEqual(true_tree,false_tree) )
    {
      print_return("",true_tree);
    } else {
      treevars_usedcnt[ATgetInt((ATermInt) ATgetArgument((ATermAppl) M,1))]++;
      print_return("",ATmakeAppl3(afunM,ATgetArgument((ATermAppl) M,0),(ATerm) true_tree,(ATerm) false_tree));
    }
  } else if ( !ATisEmpty(pars.Flist) ) {
    ATermList F = ATLgetFirst(pars.Flist);
    ATermAppl true_tree,false_tree;

    ATermList newupstack = pars.upstack;
    ATermList l = ATmakeList0();

    for (; !ATisEmpty(pars.Flist); pars.Flist=ATgetNext(pars.Flist))
    {
      if ( ATisEqual(ATgetFirst(ATLgetFirst(pars.Flist)),ATgetFirst(F)) )
      {
        newupstack = ATinsert(newupstack, (ATerm) ATgetNext(ATLgetFirst(pars.Flist)));
      } else {
        l = ATinsert(l,ATgetFirst(pars.Flist));
      }
    }

    pars.Flist = l;
    false_tree = build_tree(pars,i);
    pars.Flist = ATmakeList0();
    pars.upstack = newupstack;
    true_tree = build_tree(pars,i);

    if ( ATisEqual(true_tree,false_tree) )
    {
      print_return("",true_tree);
    } else {
      print_return("",ATmakeAppl3(afunF,ATgetArgument(ATAgetFirst(F),0),(ATerm) true_tree,(ATerm) false_tree));
    }
  } else if ( !ATisEmpty(pars.upstack) ) {
    ATermList l;

    ATermAppl r = NULL;
    ATermList readies = ATmakeList0();

    pars.stack = ATinsert(pars.stack,(ATerm) ATmakeList0());
    l = pars.upstack;
    pars.upstack = ATmakeList0();
    add_to_build_pars(&pars,l,&r,&readies);


    if ( r == NULL )
    {
      ATermAppl t = build_tree(pars,i);

      for (; !ATisEmpty(readies); readies=ATgetNext(readies))
      {
        inc_usedcnt((ATermList) ATgetArgument(ATAgetFirst(readies),2));
        inc_usedcnt((ATermList) ATgetArgument(ATAgetFirst(readies),3));
        t = ATmakeAppl3(afunC,ATgetArgument(ATAgetFirst(readies),0),(ATerm) ATmakeAppl1(afunR,ATgetArgument(ATAgetFirst(readies),1)),(ATerm) t);
      }

      print_return("",t);
    } else {
      inc_usedcnt((ATermList) ATgetArgument(r,1));
      print_return("",ATmakeAppl1(afunR,ATgetArgument(r,0)));
    }
  } else {
    if ( ATisEmpty(ATLgetFirst(pars.stack)) )
    {
      if ( ATisEmpty(ATgetNext(pars.stack)) )
      {
        print_return("",ATmakeAppl0(afunX));
      } else {
        pars.stack = ATgetNext(pars.stack);
//				print_return("",ATmakeAppl1(afunD,(ATerm) build_tree(pars,i)));
        print_return("",build_tree(pars,i));
      }
    } else {
      ATermList l = ATLgetFirst(pars.stack);
      ATermAppl r = NULL;
      ATermList readies = ATmakeList0();

      pars.stack = ATinsert(ATgetNext(pars.stack),(ATerm) ATmakeList0());
      add_to_build_pars(&pars,l,&r,&readies);

      ATermAppl tree;
      if ( r == NULL )
      {
        tree = build_tree(pars,i);
        for (; !ATisEmpty(readies); readies=ATgetNext(readies))
        {
          inc_usedcnt((ATermList) ATgetArgument(ATAgetFirst(readies),2));
          inc_usedcnt((ATermList) ATgetArgument(ATAgetFirst(readies),3));
          tree = ATmakeAppl3(afunC,ATgetArgument(ATAgetFirst(readies),0),(ATerm) ATmakeAppl1(afunR,ATgetArgument(ATAgetFirst(readies),1)),(ATerm) tree);
        }
      } else {
        inc_usedcnt((ATermList) ATgetArgument(r,1));
        tree = ATmakeAppl(afunR,ATgetArgument(r,0));
      }

      print_return("",ATmakeAppl1(afunN,(ATerm) tree));
    }
  }
}

static ATermAppl optimise_tree_aux(ATermAppl tree, ATermList stored, size_t len, size_t *max)
{
  if ( isS(tree) )
  {
    if ( len+2 > *max )
    {
      *max = len+2;
    }
    return ATmakeAppl2(afunS,ATgetArgument(tree,0),(ATerm) optimise_tree_aux(ATAgetArgument(tree,1),ATinsert(stored,ATgetArgument(tree,0)),len+1,max));
  } else if ( isM(tree) )
  {
    size_t n=len-ATindexOf(stored,ATgetArgument(tree,0),0);
    assert(n < ((size_t)1)<<(sizeof(int)*8-1));
    return ATmakeAppl3(afunM,(ATerm) ATmakeInt((int)n),(ATerm) optimise_tree_aux(ATAgetArgument(tree,1),stored,len,max),(ATerm) optimise_tree_aux(ATAgetArgument(tree,2),stored,len,max));
  } else if ( isF(tree) )
  {
    return ATmakeAppl3(afunF,ATgetArgument(tree,0),(ATerm) optimise_tree_aux(ATAgetArgument(tree,1),stored,len,max),(ATerm) optimise_tree_aux(ATAgetArgument(tree,2),stored,len,max));
  } else if ( isN(tree) )
  {
    return ATmakeAppl1(afunN,(ATerm) optimise_tree_aux(ATAgetArgument(tree,0),stored,len,max));
  } else if ( isD(tree) )
  {
    return optimise_tree_aux(ATAgetArgument(tree,0),stored,len,max);
  } else if ( isC(tree) )
  {
    return ATmakeAppl3(afunC,ATgetArgument(tree,0),(ATerm) optimise_tree_aux(ATAgetArgument(tree,1),stored,len,max),(ATerm) optimise_tree_aux(ATAgetArgument(tree,2),stored,len,max));
  } else if ( isR(tree) )
  {
    return tree;
  } else /* isX(tree) */ {
    assert(isX(tree));
    return tree;
  }
}

static ATermAppl optimise_tree(ATermAppl tree,size_t *max)
{
  return optimise_tree_aux(tree,ATmakeList0(),ATERM_NON_EXISTING_POSITION,max);
}

#ifdef _INNER_STORE_TREES
ATermAppl RewriterInnermost::create_tree(ATermList rules, size_t opid, size_t *max_vars, ATermInt trueint)
#else
static ATermAppl create_tree(ATermList rules, size_t /*opid*/, size_t *max_vars, ATermInt trueint)
#endif
  // Create a match tree for OpId int2term[opid] and update the value of
  // *max_vars accordingly.
  //
  // Pre:  rules is a list of rewrite rules for int2term[opid] in the
  //       INNER internal format
  //       opid is a valid entry in int2term
  //       max_vars is a valid pointer to an integer
  // Post: *max_vars is the maximum of the original *max_vars value and
  //       the number of variables in the result tree
  // Ret:  A match tree for int2term[opid]
{
//gsfprintf(stderr,"%P (%i)\n",int2term[opid],opid);
  // Create sequences representing the trees for each rewrite rule and
  // store the total number of variables used in these sequences.
  // (The total number of variables in all sequences should be an upper
  // bound for the number of variable in the final tree.)
  ATermList rule_seqs = ATmakeList0();
  size_t total_rule_vars = 0;
  for (; !ATisEmpty(rules); rules=ATgetNext(rules))
  {
    rule_seqs = ATinsert(rule_seqs, (ATerm) create_sequence((ATermAppl) ATgetFirst(rules),&total_rule_vars, trueint));
  }

  // Generate initial parameters for built_tree
  build_pars init_pars;
  ATermAppl r = NULL;
  ATermList readies = ATmakeList0();

  initialise_build_pars(&init_pars);
  add_to_build_pars(&init_pars,rule_seqs,&r,&readies);

  ATermAppl tree;
  if ( r == NULL )
  {
    MCRL2_SYSTEM_SPECIFIC_ALLOCA(a,size_t,total_rule_vars);
    treevars_usedcnt = a;
//		treevars_usedcnt = (size_t *) malloc(total_rule_vars*sizeof(size_t));
    tree = build_tree(init_pars,0);
//		free(treevars_usedcnt);
    for (; !ATisEmpty(readies); readies=ATgetNext(readies))
    {
      tree = ATmakeAppl3(afunC,ATgetArgument(ATAgetFirst(readies),0),(ATerm) ATmakeAppl1(afunR,ATgetArgument(ATAgetFirst(readies),1)),(ATerm) tree);
    }
  } else {
    tree = ATmakeAppl1(afunR,ATgetArgument(r,0));
  }
  //ATfprintf(stderr,"tree(%t): %t\n",int2term[opid],tree);

  finalise_build_pars(&init_pars);

  size_t max_tree_vars = 0;
  tree = optimise_tree(tree,&max_tree_vars);
  if ( max_tree_vars > *max_vars )
  {
    *max_vars = max_tree_vars;
  }

#ifdef _INNER_STORE_TREES
  char s[100],t[100];
  sprintf(s,"tree_%i_%s",opid,ATgetName(ATgetAFun(ATAgetArgument(int2term[opid],0))));
  sprintf(t,"tree_%i_%s.dot",opid,ATgetName(ATgetAFun(ATAgetArgument(int2term[opid],0))));
  tree2dot(tree,s,t);
#endif

  return tree;
}

//#define TMA_DEBUG
ATermList RewriterInnermost::tree_matcher_aux(ATerm t, ATermAppl *tree, ATermAppl *vars, ATerm *vals, size_t *len)
{
  ATermList args;

  args = ATmakeList0();
  if ( ATisList(t) )
  {
    args = ATgetNext((ATermList) t);
#ifdef TMA_DEBUG
    t = ATgetFirst((ATermList) t);
#endif
  }

#ifdef TMA_DEBUG
ATfprintf(stderr,"aux: %t %t\n",t,args);
#endif
  while ( !ATisEmpty(args) || isD(*tree) || isC(*tree) )
  {
#ifdef TMA_DEBUG
ATfprintf(stderr,"loop: %t\n",args);
#endif
    if ( isN(*tree) )
    {
#ifdef TMA_DEBUG
ATfprintf(stderr,"N\n");
#endif
      args = ATgetNext(args);
      if ( !ATisEmpty(args) )
      {
        *tree = ATAgetArgument(*tree,0);
      }
    } else if ( isS(*tree) )
    {
#ifdef TMA_DEBUG
ATfprintf(stderr,"S %t\n",ATgetArgument(*tree,0));
#endif
      vars[*len] = (ATermAppl) ATgetArgument(*tree,0);
      vals[*len] = ATgetFirst(args);
      (*len)++;
      *tree = ATAgetArgument(*tree,1);
    } else if ( isR(*tree) )
    {
#ifdef TMA_DEBUG
ATfprintf(stderr,"R\n");
#endif
      return args;
    } else if ( isM(*tree) )
    {
#ifdef TMA_DEBUG
ATfprintf(stderr,"M %t\n",ATgetArgument(*tree,0));
#endif
      if ( ATisEqual(ATgetFirst(args),vals[ATgetInt((ATermInt) ATgetArgument(*tree,0))]) )
      {
        *tree = ATAgetArgument(*tree,1);
      } else {
        *tree = ATAgetArgument(*tree,2);
      }
/*		} else if ( isD(*tree) )
    {
#ifdef TMA_DEBUG
ATfprintf(stderr,"D\n");
#endif
      *tree = ATAgetArgument(*tree,0);
#ifdef TMA_DEBUG
ATfprintf(stderr,"return\n");
#endif
      return args;*/
    } else if ( isF(*tree) )
    {
#ifdef TMA_DEBUG
ATfprintf(stderr,"F %t\n",ATgetArgument(*tree,0));
#endif
      if ( ATisEqual(ATgetFirst(args),ATgetArgument(*tree,0)) || (ATisList(ATgetFirst(args)) && ATisEqual(ATgetFirst(ATLgetFirst(args)),ATgetArgument(*tree,0))) )
      {
#ifdef TMA_DEBUG
ATfprintf(stderr,"true\n");
#endif
        *tree = ATAgetArgument(*tree,1);
        tree_matcher_aux(ATgetFirst(args),tree,vars,vals,len);
      } else {
#ifdef TMA_DEBUG
ATfprintf(stderr,"false\n");
#endif
        *tree = ATAgetArgument(*tree,2);
      }
    } else if ( isX(*tree) )
    {
#ifdef TMA_DEBUG
ATfprintf(stderr,"X\n");
#endif
      return args;
    } else if ( isC(*tree) )
    {
#ifdef TMA_DEBUG
ATfprintf(stderr,"C\n");
#endif
      if ( ATisEqual(build(ATgetArgument(*tree,0),ATERM_NON_EXISTING_POSITION,vars,vals,*len),trueint) )
      {
#ifdef TMA_DEBUG
ATfprintf(stderr,"true (c)\n");
#endif
        *tree = ATAgetArgument(*tree,1);
      } else {
#ifdef TMA_DEBUG
ATfprintf(stderr,"false (c)\n");
#endif
        *tree = ATAgetArgument(*tree,2);
      }
    } else {
#ifdef TMA_DEBUG
ATfprintf(stderr,"? %t\n",*tree);
#endif
      assert(false);
    }
  }
#ifdef TMA_DEBUG
ATfprintf(stderr,"no more args\n");
#endif
  return args;
}

ATerm RewriterInnermost::tree_matcher(ATermList t, ATermAppl tree)
{
  MCRL2_SYSTEM_SPECIFIC_ALLOCA(vars,ATermAppl,max_vars);
  MCRL2_SYSTEM_SPECIFIC_ALLOCA(vals,ATerm,max_vars);
  size_t len = 0;

  while ( isC(tree) )
  {
    if ( ATisEqual(build(ATgetArgument(tree,0),ATERM_NON_EXISTING_POSITION,vars,vals,len),trueint) )
    {
      tree = (ATermAppl) ATgetArgument(tree,1); // Was 0????  JFG This was a very odd error.
    } else {
      tree = (ATermAppl) ATgetArgument(tree,2); // Was 1????  JFG This was a very odd error.
    }
  }

  ATermList rargs;
  if ( isR(tree) )
  {
    rargs = ATgetNext((ATermList) t);
  } else {
    rargs = tree_matcher_aux((ATerm) t,&tree,vars,vals,&len);
    rargs = ATgetNext(rargs);
  }

  if ( isR(tree) )
  {

    ATerm rslt = ATgetArgument(tree,0);
    size_t rslt_len;
    if ( ATisList(rslt) )
    {
      rslt_len = ATgetLength((ATermList) rslt)-1;
      if ( !ATisEmpty(rargs) )
      {
        rslt = (ATerm) ATconcat((ATermList) rslt,rargs);
      }
    } else {
      rslt_len = 0;
      if ( !ATisEmpty(rargs) )
      {
        rslt = (ATerm) ATinsert(rargs,rslt);
      }
    }

    ATerm r = build(rslt,rslt_len,vars,vals,len);

    return r;
  } else {
    return NULL;
  }
}

ATerm RewriterInnermost::OpId2Int(ATermAppl Term, bool add_opids)
{
  ATermInt i;

  if ( (i = (ATermInt) ATtableGet(term2int,(ATerm) Term)) == NULL )
  {
    if ( !add_opids )
    {
      return (ATerm) Term;
    }
    assert(num_opids<((size_t)1)<<(sizeof(int)*8-1));
    i = ATmakeInt((int)num_opids);
    ATtablePut(term2int,(ATerm) Term,(ATerm) i);
    num_opids++;
  }

  return (ATerm) i;
}

ATerm RewriterInnermost::toInner(ATermAppl Term, bool add_opids)
{
  ATermList l;

  if ( !gsIsDataAppl(Term) )
  {
    if ( gsIsOpId(Term) )
    {
      return (ATerm) OpId2Int(Term,add_opids);
    } else {
      return (ATerm) Term;
    }
  }

  l = ATmakeList0();
        if ( gsIsDataAppl(Term) )
        {
                for ( ATermList args = ATLgetArgument((ATermAppl) Term,1) ; !ATisEmpty(args) ; args = ATgetNext(args))
                {
                        l = ATinsert(l,(ATerm) toInner((ATermAppl) ATgetFirst(args),add_opids));
                }

                l = ATreverse(l);

                ATerm arg0 = toInner(ATAgetArgument((ATermAppl) Term,0),add_opids);

                if ( ATisList(arg0) )
                {
                        l = ATconcat((ATermList) arg0, (ATermList) l);
                }
                else
                {
                        l = ATinsert(l, (ATerm) arg0);
                }
        }

  return (ATerm) l;
}

ATermAppl RewriterInnermost::fromInner(ATerm Term)
{
        assert(Term!=NULL);

  if ( !ATisList(Term) )
  {
    if ( ATisInt(Term) )
    {
      return int2term[ATgetInt((ATermInt) Term)];
    } else {
      return (ATermAppl) Term;
    }
  }

  if ( ATisEmpty((ATermList) Term) )
  {
    throw mcrl2::runtime_error(std::string(NAME) + ": invalid inner format term.");
  }

        //Reconstruct term structure
  ATermList l = (ATermList) Term;
  ATerm t = ATgetFirst(l);
  ATermAppl a;
  if ( ATisInt(t) )
  {
    a = int2term[ATgetInt((ATermInt) t)];
  } else {
    a = (ATermAppl) t;
  }


  l = ATgetNext(l);
        if(gsIsOpId(a) || gsIsDataVarId(a))
        {
                ATermAppl sort = ATAgetArgument(a, 1);
                while(is_function_sort(sort_expression(sort)) && !ATisEmpty(l))
                {
                        ATermList sort_dom = ATLgetArgument(sort, 0);
                        ATermList list = ATmakeList0();
                        while (!ATisEmpty(sort_dom))
                        {
                                list = ATinsert(list, (ATerm) fromInner(ATgetFirst(l)));
                                sort_dom = ATgetNext(sort_dom);
                                l = ATgetNext(l);
                        }
                        list = ATreverse(list);
                        a = gsMakeDataAppl(a, list);
                        sort = ATAgetArgument(sort, 1);
                }
        }

  return a;
}

RewriterInnermost::RewriterInnermost(const data_specification &DataSpec)
{
  ATermList m,n;
  ATermTable tmp_eqns;
  ATermInt i;

  initialise_common();

  num_opids = 0;
  max_vars = 0;

  tmp_eqns = ATtableCreate(100,50); // XXX would be nice to know the number op OpIds
  term2int = ATtableCreate(100,50);

  trueint = NULL;
  ATprotectInt(&trueint);
  trueint = (ATermInt) OpId2Int(sort_bool::true_(),true);

  /*l = opid_eqns;
  for (; !ATisEmpty(l); l=ATgetNext(l))
  {
    // XXX only adds the last rule where lhs is an opid; this might go "wrong" if this rule is removed later
//		ATtablePut(tmp_eqns,OpId2Int(ATAgetArgument(ATAgetFirst(l),2),true),(ATerm) ATmakeList1((ATerm) ATmakeList4((ATerm) ATmakeList0(),toInner(ATAgetArgument(ATAgetFirst(l),1),true),(ATerm) ATmakeList0(),toInner(ATAgetArgument(ATAgetFirst(l),3),true))));
    ATtablePut(tmp_eqns,OpId2Int(ATAgetArgument(ATAgetFirst(l),2),true),(ATerm) ATmakeList1((ATerm) ATmakeAppl4(ruleAFun,(ATerm) ATmakeList0(),toInner(ATAgetArgument(ATAgetFirst(l),1),true),(ATerm) ATmakeList0(),toInner(ATAgetArgument(ATAgetFirst(l),3),true))));
  }*/

//	l = dataappl_eqns;
  const data_specification::equations_const_range l = DataSpec.equations();
  for (atermpp::set< data_equation >::const_iterator j=l.begin(); j!=l.end(); ++j)
  {
    try
    {
      CheckRewriteRule(*j);
    } catch ( std::runtime_error &e ) {
      gsWarningMsg("%s\n",e.what());
      continue;
    }

    ATerm u = toInner(j->lhs(),true);
    ATerm head;
    ATermList args;

    if ( ATisInt(u) )
    {
      head = u;
      args = ATmakeList0();
    } else {
      head = ATgetFirst((ATermList) u);
      args = ATgetNext((ATermList) u);
    }

    if ( (n = (ATermList) ATtableGet(tmp_eqns,head)) == NULL )
    {
      n = ATmakeList0();
    }
    n = ATinsert(n,(ATerm) ATmakeAppl4(ruleAFun,
                                       (ATerm) static_cast<ATermList>(j->variables()),
                                       (ATerm) toInner(j->condition(),true),
                                       (ATerm) args,
                                       (ATerm) toInner(j->rhs(),true)));
    ATtablePut(tmp_eqns,head,(ATerm) n);
  }

  int2term = (ATermAppl *) malloc(num_opids*sizeof(ATermAppl));
  inner_eqns = (ATermList *) malloc(num_opids*sizeof(ATermList));
  inner_trees = (ATermAppl *) malloc(num_opids*sizeof(ATermAppl));
  for (size_t i=0; i < num_opids; i++)
  {
    int2term[i] = NULL;
    inner_eqns[i] = NULL;
    inner_trees[i] = NULL;
  }
  ATprotectArray((ATerm *) int2term,num_opids);
  ATprotectArray((ATerm *) inner_eqns,num_opids);
  ATprotectArray((ATerm *) inner_trees,num_opids);

  ATermList l1 = ATtableKeys(term2int);
  for (; !ATisEmpty(l1); l1=ATgetNext(l1))
  {
    i = (ATermInt) ATtableGet(term2int,ATgetFirst(l1));
    int2term[ATgetInt(i)] = ATAgetFirst(l1);
  }

  l1 = ATtableKeys(term2int);
  for (; !ATisEmpty(l1); l1=ATgetNext(l1))
  {
    i = (ATermInt) ATtableGet(term2int,ATgetFirst(l1));
    if ( (m = (ATermList) ATtableGet(tmp_eqns,(ATerm) i)) != NULL )
    {
      inner_eqns[ATgetInt(i)] = m;
      inner_trees[ATgetInt(i)] = create_tree(m,ATgetInt(i),&max_vars, trueint);
    }
  }

  ATtableDestroy(tmp_eqns);

  need_rebuild = false;
}

RewriterInnermost::~RewriterInnermost()
{
  ATtableDestroy(term2int);
  ATunprotectInt(&trueint);
  ATunprotectArray((ATerm *) int2term);
  ATunprotectArray((ATerm *) inner_eqns);
  ATunprotectArray((ATerm *) inner_trees);
  free(int2term);
  free(inner_eqns);
  free(inner_trees);

  finalise_common();
}

bool RewriterInnermost::addRewriteRule(ATermAppl Rule)
{
  ATermAppl m;
  ATermInt j;

  try
  {
    CheckRewriteRule(Rule);
  } catch ( std::runtime_error &e ) {
    gsWarningMsg("%s\n",e.what());
    return false;
  }

  size_t old_num = num_opids;

  ATermAppl a = ATAgetArgument(Rule,2);
  if ( gsIsOpId(a) )
  {
    j = (ATermInt) OpId2Int(a,true);
    m = ATmakeAppl4(ruleAFun,(ATerm) ATmakeList0(),toInner(ATAgetArgument(Rule,1),true),(ATerm) ATmakeList0(),toInner(ATAgetArgument(Rule,3),true));
  } else {
    ATermList l = (ATermList) toInner(a,true);
    j = (ATermInt) ATgetFirst(l);
    m = ATmakeAppl4(ruleAFun,ATgetArgument(Rule,0),toInner(ATAgetArgument(Rule,1),true),(ATerm) ATgetNext(l),toInner(ATAgetArgument(Rule,3),true));
  }

  if ( num_opids > old_num )
  {
    ATunprotectArray((ATerm *) int2term);
    ATunprotectArray((ATerm *) inner_eqns);
    ATunprotectArray((ATerm *) inner_trees);

    int2term = (ATermAppl *) realloc(int2term,num_opids*sizeof(ATermAppl));
    inner_eqns = (ATermList *) realloc(inner_eqns,num_opids*sizeof(ATermList));
    inner_trees = (ATermAppl *) realloc(inner_trees,num_opids*sizeof(ATermAppl));
    for (size_t k=old_num; k < num_opids; k++)
    {
      int2term[k] = NULL;
      inner_eqns[k] = NULL;
      inner_trees[k] = NULL;
    }
    ATprotectArray((ATerm *) int2term,num_opids);
    ATprotectArray((ATerm *) inner_eqns,num_opids);
    ATprotectArray((ATerm *) inner_trees,num_opids);

    ATermList l = ATtableKeys(term2int);
    for (; !ATisEmpty(l); l=ATgetNext(l))
    {
      ATermInt i = (ATermInt) ATtableGet(term2int,ATgetFirst(l));
      if ( (size_t)ATgetInt(i) >= old_num )
      {
        int2term[ATgetInt(i)] = ATAgetFirst(l);
      }
    }
  }

  if ( inner_eqns[ATgetInt(j)] == NULL )
  {
    inner_eqns[ATgetInt(j)] = ATmakeList1((ATerm) m);
  } else {
    inner_eqns[ATgetInt(j)] = ATinsert(inner_eqns[ATgetInt(j)],(ATerm) m);
  }

  need_rebuild = true;

  return true;
}

bool RewriterInnermost::removeRewriteRule(ATermAppl Rule)
{
  ATermAppl m;
  ATerm t;

  ATermAppl a = ATAgetArgument(Rule,2);
  if ( gsIsOpId(a) )
  {
    t = OpId2Int(a,false);
    m = ATmakeAppl4(ruleAFun,(ATerm) ATmakeList0(),toInner(ATAgetArgument(Rule,1),true),(ATerm) ATmakeList0(),toInner(ATAgetArgument(Rule,3),true));
  } else {
    ATermList l = (ATermList) toInner(a,false);
    t = ATgetFirst(l);
    m = ATmakeAppl4(ruleAFun,ATgetArgument(Rule,0),toInner(ATAgetArgument(Rule,1),true),(ATerm) ATgetNext(l),toInner(ATAgetArgument(Rule,3),true));
  }

  if ( ATisInt(t) )
  {
    ATermList l = inner_eqns[ATgetInt((ATermInt) t)];
    ATermList n = ATmakeList0();
    for (; !ATisEmpty(l); l=ATgetNext(l))
    {
      if ( !ATisEqual(m,ATgetFirst(l)) )
      {
        ATinsert(n,ATgetFirst(l));
      }
    }
    if ( ATisEmpty(n) )
    {
      inner_eqns[ATgetInt((ATermInt) t)] = NULL;
    } else {
      inner_eqns[ATgetInt((ATermInt) t)] = n;
    }

    need_rebuild = true;
  }

  return true;
}

ATermList RewriterInnermost::build_args(ATermList args, size_t buildargs, ATermAppl *vars, ATerm *vals, size_t len)
{
  if ( (buildargs == 0) || ATisEmpty(args) )
  {
    return args;
  } else {
    return ATinsert(build_args(ATgetNext(args),buildargs-1,vars,vals,len),build(ATgetFirst(args),ATERM_NON_EXISTING_POSITION,vars,vals,len));
  }
}

ATerm RewriterInnermost::build(ATerm Term, size_t buildargs, ATermAppl *vars, ATerm *vals, size_t len)
{
  if ( ATisList(Term) )
  {
    ATerm head = ATgetFirst((ATermList) Term);
    ATermList args = ATgetNext((ATermList) Term);

    if ( buildargs == ATERM_NON_EXISTING_POSITION )
    {
      buildargs = ATgetLength(args);
    }

    args = build_args(args,buildargs,vars,vals,len);

    size_t b = 1;
    while ( !ATisInt(head) && b )
    {
      for (size_t i=0; i<len; i++)
      {
        if ( ATisEqual(head,vars[i]) )
        {
          if ( ATisList(vals[i]) )
          {
            head = ATgetFirst((ATermList) vals[i]);
            args = ATconcat(ATgetNext((ATermList) vals[i]),args);
          } else {
            head = vals[i];
          }
          break;
        }
        b = 0;
      }
    }
    if ( ATisInt(head) )
    {
      return rewrite_func((ATermInt) head,args);
    } else {
      return (ATerm) ATinsert(args,head);
    }
  } else if ( ATisInt(Term) ) {
    return rewrite_func((ATermInt) Term, ATmakeList0());
  } else {
    for (size_t i=0; i<len; i++)
    {
      if ( ATisEqual(Term,vars[i]) )
      {
        return vals[i];
      }
    }
    return Term;
  }
}

ATerm RewriterInnermost::rewrite_func(ATermInt op, ATermList args)
{
  ATermAppl tree;
//gsfprintf(stderr,"rewrite_func(%T,%T)\n\n",op,args);
//gsfprintf(stderr,"rewrite_func(%P)\n\n",fromInner((ATerm) ATinsert(args,(ATerm) op)));

  if ( (tree = inner_trees[ATgetInt(op)]) != NULL )
  {
    ATerm r = tree_matcher(ATinsert(args,(ATerm) op),tree);

    if ( r != NULL )
    {
      return r;
    }
  }

  if ( ATisEmpty(args) )
  {
    return (ATerm) op;
  } else {
    return (ATerm) ATinsert(args,(ATerm) op);
  }
}

ATermList RewriterInnermost::rewrite_listelts(ATermList l)
{
  if ( ATisEmpty(l) )
  {
    return l;
  } else {
    return ATinsert(rewrite_listelts(ATgetNext(l)),rewrite_aux(ATgetFirst(l)));
  }
}

ATerm RewriterInnermost::rewrite_aux(ATerm Term)
{
//gsfprintf(stderr,"rewrite_aux(%T)\n\n",Term);
//gsfprintf(stderr,"rewrite_aux(%P)\n\n",fromInner(Term));
  if ( ATisList(Term) )
  {
    ATerm head = ATgetFirst((ATermList) Term);
    ATermList l = ATgetNext((ATermList) Term);

    // rewrite arguments
    l = rewrite_listelts(l);

    // head is should be a size_t or a var
    // if it's a var, see if it needs to be substituted
    if ( !ATisInt(head) )
    {
      if ( !gsIsOpId(head) ) // this internal rewrite format allows opids
                       // XXX change to main internal rewrite format!
      {
        assert(ATisAppl(head) && gsIsDataVarId((ATermAppl) head));

        ATerm a = lookupSubstitution((ATermAppl) head);

        // if a is a list, concatenate the arguments
        // always set head again, as 'a' might have been changed
        if ( ATisList(a) )
        {
          head = ATgetFirst((ATermList) a);
          l = ATconcat(ATgetNext((ATermList) a),l);
        } else {
          head = a;
        }
      }
    }

    // head is should be a size_t or a var
    if ( ATisInt(head) )
    {
      // head is a size_t, thus we can try to rewrite further
      Term = rewrite_func((ATermInt) head, l);
    } else {
      // head is a var, so there is nothing left to do but to
      // reattach the head
      Term = (ATerm) ATinsert(l,head);
    }

    return Term;
  } else if ( ATisInt(Term) )
  {
    return rewrite_func((ATermInt) Term, ATmakeList0());
  } else {
    return lookupSubstitution((ATermAppl) Term);
  }
}

ATerm RewriterInnermost::rewriteInternal(ATerm Term)
{
  if ( need_rebuild )
  {
    for (size_t i=0; i < num_opids; i++)
    {
      if ( (inner_trees[i] == NULL) && (inner_eqns[i] != NULL) )
      {
        inner_trees[i] = create_tree(inner_eqns[i],i,&max_vars, trueint);
      }
    }
  }

  return rewrite_aux(Term);
}

ATermAppl RewriterInnermost::rewrite(ATermAppl Term)
{
        /*gsDebugMsg("Rewriting term: %T\n\n", Term);
        gsDebugMsg("toInner(Term): %T\n\n", toInner(Term,false));
        gsDebugMsg("fromInner(toInner(Term)): %T\n\n\n", fromInner(toInner(Term,false)));*/
  return fromInner(rewriteInternal(toInner(Term,false)));
}

ATerm RewriterInnermost::toRewriteFormat(ATermAppl Term)
{
  return toInner(Term,false);
}

ATermAppl RewriterInnermost::fromRewriteFormat(ATerm Term)
{
  return fromInner(Term);
}

RewriteStrategy RewriterInnermost::getStrategy()
{
  return GS_REWR_INNER;
}

    }
  }
}

