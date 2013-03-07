#ifndef WEALTH_TEMPLATE_H
#define WEALTH_TEMPLATE_H

#include "wealth.h"
#include "line_search.Template.h"
#include <cassert>

template <class Bidder>
void
DualWealthArray::initialize_wealth_bid_array(Bidder const& f, int nRounds)
{
  std::vector< double > wealths, bids;
  int          i        (0);
  double       minW     (f.total_wealth());   // minW tracks wealth as we subtract away bids
  int          checkPt  (10000);
  while (mW0 <= minW)                         // spend down wealth until reach initial wealth
  { double b (f(i));
    wealths.push_back(minW);
    bids.push_back(b);
    minW -= b;
    ++i;
    if(i == checkPt)
    { std::clog << "WLTH: Warning. Spending down slowly. At " << i << " with minW " << minW << " > W0 =" << mW0 << std::endl;
      checkPt += 50000;
    }
  }
  for (double j=0; j<=nRounds; ++j)           // spend down from W0 by number of rounds as if lose each, with an added extra round
  { double b = f(i+j);
    wealths.push_back(minW);
    bids.push_back(b);
    minW -= b;
  }
  i = (int) wealths.size();
  while ((mW0-grid_delta(mW0)) < minW)       // spend enough to get below W0 (univ is slow)
  { double b = f(i++);
    wealths.push_back(minW);
    bids.push_back(b);
    minW -= b;
  }
  std::clog << "WLTH: Note that mW0-grid_delta = " << mW0-grid_delta(mW0) << " > " << minW << " = min wealth" <<std::endl;
  std::clog << "WLTH: Wealth vector uses " << wealths.size() << " elements, with min wealth at position "
	    << wealths.size()-1 << "=" << minW << " with bid=" << bids[bids.size()-1] << std::endl;
  // initialize wealth array
  double delta   (grid_delta(f.total_wealth()));        // spacing
  double wealth  (0.1 * floor(10*f.total_wealth()));    // round to 0.1 to hit typical W0 value on grid
  int    iStart  (0);         
  int    limit   (1000);      
  while ((minW < wealth) && limit--)  // dont step below smallest tabulated wealth
  { double indx (Line_Search::invert_monotone(wealth, iStart, wealths));
    iStart = floor(indx);
    double shr (indx-iStart);
    double bid (bids[iStart]*(1.0-shr) + bids[iStart+1]*shr);
    mWealthBid.push_back( std::make_pair(wealth,bid) );
    wealth -= delta;
    delta = grid_delta(wealth);
  }
  if (limit <= 0) std::clog << "WLTH: *** ERROR *** Ran out of space for internal array.\n" << std::endl;
  // find position with wealth equal to W0 (up to rounding)
  i = 0;
  while ((0.000001 < (mWealthBid[i].first - mW0)) && (i < (int)mWealthBid.size())) ++i;
  mZeroIndex = i;
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
