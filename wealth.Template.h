#ifndef WEALTH_TEMPLATE_H
#define WEALTH_TEMPLATE_H

#include "wealth.h"
#include "line_search.Template.h"
#include <cassert>

template <class Rule>
void
DualWealthArray::initialize_wealth_bid_array(Rule const& f, double maxWealth, int nRounds)
{
  // compute smallest attainable wealth, verify it is less than initial wealth
  double minWealth (mW0);
  for(int i=0; i<nRounds; ++i)
    minWealth -= f(minWealth);
  if (mW0 < minWealth + grid_delta(mW0))
  { minWealth = mW0 - grid_delta(mW0);
    std::clog << "WLTH: Slow decaying wealth function did not reach enough below mW0; min set to " << minWealth << std::endl;
  }
  // initialize wealth array
  double     delta   (grid_delta(maxWealth));               // tracks spacing
  double     wealth  (0.1 * floor(10*maxWealth));           // round to 0.1 to hit typical W0 value on grid
  int        limit   (1500);
  do
  { double bid (f(wealth));
    mWealthBid.push_back( std::make_pair(wealth,bid) );
    wealth -= delta;
    delta = grid_delta(wealth);
  } while ((minWealth <= wealth+0.000001) && limit--);               
  if (limit <= 0)
  { std::clog << "WLTH: *** ERROR *** Ran out of space for internal array.\n" << std::endl;
    assert(false);
  }
  // find mZeroIndex position for which wealth equals W0 (up to rounding)
  { int i (0);
    while ((0.00001 < (mWealthBid[i].first - mW0)) && (i < (int)mWealthBid.size())) ++i;
    mZeroIndex = i;
  }
  if (((int)mWealthBid.size()-1) == mZeroIndex)
  { std::clog << "WLTH: *** Error *** W0 position in the last element of wealth grid (size=" << mWealthBid.size()<<") at " << mZeroIndex << std::endl;
    assert (false);
  }
  for(int i=0; i < (int)mWealthBid.size(); ++i)
    if (1.0 < mWealthBid[i].second)
    { std::clog << "WLTH: Found dual wealth bid " << mWealthBid[i].second << " @ " << i << " larger than 1; reduced to 1." << std::endl;
      mWealthBid[i].second = 1.0;
    }
}

#endif
