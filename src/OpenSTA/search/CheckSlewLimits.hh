// OpenSTA, Static Timing Analyzer
// Copyright (c) 2020, Parallax Software, Inc.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#pragma once

#include "MinMax.hh"
#include "Transition.hh"
#include "NetworkClass.hh"
#include "GraphClass.hh"
#include "Delay.hh"
#include "SdcClass.hh"

namespace sta {

class StaState;
class DcalcAnalysisPt;
class Corner;

class CheckSlewLimits
{
public:
  CheckSlewLimits(const StaState *sta);
  // corner=nullptr checks all corners.
  void checkSlew(const Pin *pin,
		 const Corner *corner,
		 const MinMax *min_max,
		 bool check_clks,
		 // Return values.
		 // Corner is nullptr for no slew limit.
		 const Corner *&corner1,
		 const RiseFall *&rf,
		 Slew &slew,
		 float &limit,
		 float &slack) const;
  // corner=nullptr checks all corners.
  PinSeq *pinSlewLimitViolations(const Corner *corner,
				 const MinMax *min_max);
  // corner=nullptr checks all corners.
  Pin *pinMinSlewLimitSlack(const Corner *corner,
			    const MinMax *min_max);

protected:
  void checkSlews1(const Pin *pin,
		   const Corner *corner,
		   const MinMax *min_max,
		   bool check_clks,
		   // Return values.
		   const Corner *&corner1,
		   const RiseFall *&rf,
		   Slew &slew,
		   float &limit,
		   float &slack) const;
  void checkSlews1(Vertex *vertex,
		   const Corner *corner1,
		   const MinMax *min_max,
		   bool check_clks,
		   // Return values.
		   const Corner *&corner,
		   const RiseFall *&rf,
		   Slew &slew,
		   float &limit,
		   float &slack) const;
  void checkSlew(Vertex *vertex,
		 const Corner *corner1,
		 const RiseFall *rf1,
		 const MinMax *min_max,
		 float limit1,
		 // Return values.
		 const Corner *&corner,
		 const RiseFall *&rf,
		 Slew &slew,
		 float &slack,
		 float &limit) const;
  void findLimit(const Pin *pin,
		 const Vertex *vertex,
		 const RiseFall *rf,
		 const MinMax *min_max,
		 bool check_clks,
		 // Return values.
		 float &limit1,
		 bool &limit1_exists) const;
  void pinSlewLimitViolations(Instance *inst,
			      const Corner *corner,
			      const MinMax *min_max,
			      PinSeq *violators);
  void pinMinSlewLimitSlack(Instance *inst,
			    const Corner *corner,
			    const MinMax *min_max,
			    // Return values.
			    Pin *&min_slack_pin,
			    float &min_slack);
  void clockDomains(const Vertex *vertex,
		    // Return value.
		    ClockSet &clks) const;

  const StaState *sta_;
};

} // namespace
