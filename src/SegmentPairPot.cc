// Copyright 2008, 2010 Martin C. Frith

#include "SegmentPairPot.hh"
#include <cassert>

namespace cbrc{

void SegmentPairPot::sort(){
  std::sort( items.begin(), items.end(), itemLess );

  // Remove duplicates.  We assume that, after sorting, duplicates are
  // consecutive.  This will be true if non-duplicates never overlap,
  // which is true if all the SegmentPairs are "optimal".
  iterator newEnd = std::unique( items.begin(), items.end() );
  items.erase( newEnd, items.end() );

  iters.clear();

  for( iterator i = items.begin(); i < items.end(); ++i ){
    iters.push_back(i);
  }

  std::sort( iters.begin(), iters.end(), iterLess );
}

void SegmentPairPot::markOverlaps( const SegmentPair& sp ){
  iterator i = std::lower_bound( items.begin(), items.end(), sp, itemLess );

  // Assume we need to check just one previous item.  This assumption
  // will be true provided that the items never overlap each other.
  if( i > items.begin() &&
      (i-1)->diagonal() == sp.diagonal() &&
      (i-1)->end1() > sp.beg1() ){
    mark( *(i-1) );
  }

  while( i < items.end() &&
	 i->diagonal() == sp.diagonal() &&
	 i->beg1() < sp.end1() ){
    mark( *i );
    ++i;
  }
}

void SegmentPairPot::markAllOverlaps( const std::vector<SegmentPair>& sps ){
  for( const_iterator i = sps.begin(); i < sps.end(); ++i ){
    markOverlaps( *i );
  }
}

void SegmentPairPot::markTandemRepeats( const SegmentPair& sp,
					indexT maxDistance ){
  assert( !items.empty() );

  // Careful: if we are self-comparing lots of short sequences, there
  // may be many items on the same diagonal as sp.

  SegmentPair nextDiagonal( sp.beg1() + 1, sp.beg2(), 0, 0 );
  iterator n = std::lower_bound( items.begin(), items.end(),
				 nextDiagonal, itemLess );
  if( n == items.end() )  n = items.begin();

  iterator j = n;
  do{  // funny loop to deal with wrap-around
    if( j->diagonal() - sp.diagonal() > maxDistance )  break;
    if( j->beg2() >= sp.beg2() && j->end1() <= sp.end1() )  mark( *j );
    ++j;
    if( j == items.end() )  j = items.begin();
  }while( j != n );

  SegmentPair prevDiagonal( sp.end1() - 1, sp.end2(), 0, 0 );
  iterator p = std::lower_bound( items.begin(), items.end(),
				 prevDiagonal, itemLess );
  if( p == items.end() )  p = items.begin();

  iterator k = p;
  do{  // funny loop to deal with wrap-around
    if( k == items.begin() )  k = items.end();
    --k;
    if( sp.diagonal() - k->diagonal() > maxDistance )  break;
    if( k->beg1() >= sp.beg1() && k->end2() <= sp.end2() )  mark( *k );
  }while( k != p );
}

}
