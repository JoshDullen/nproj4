//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Stub projection implementation for lat/long coordinates. We
//			don't actually change the coordinates, but we want proj=latlong
//			to act sort of like a projection.
// Author:	Frank Warmerdam, warmerda@home.com
//
//*****************************************************************************
// Copyright (c) 2000, Frank Warmerdam
// Copyright (c) 2008-2011 by the Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//*****************************************************************************

// very loosely based upon DMA code by Bradford W. Drew
namespace Free.Ports.Proj4.Projections
{
	class PJ_lonlat : PJ
	{
		public override string Name { get { return "lonlat"; } }
		public override string DescriptionName { get { return "Lat/long (Geodetic)"; } }
		public override string DescriptionType { get { return ""; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		protected XY forward(LP lp)
		{
			XY xy;
			xy.x=lp.lam/a;
			xy.y=lp.phi/a;
			return xy;
		}

		protected LP inverse(XY xy)
		{
			LP lp;
			lp.phi=xy.y*a;
			lp.lam=xy.x*a;
			return lp;
		}

		public override PJ Init()
		{
			is_latlong=true;
			x0=0.0;
			y0=0.0;
			inv=inverse;
			fwd=forward;

			return this;
		}
	}

	class PJ_latlon : PJ_lonlat
	{
		public override string Name { get { return "latlon"; } }
		public override string DescriptionName { get { return "Lat/long (Geodetic alias)"; } }
	}

	class PJ_latlong : PJ_lonlat
	{
		public override string Name { get { return "latlong"; } }
		public override string DescriptionName { get { return "Lat/long (Geodetic alias)"; } }
	}

	class PJ_longlat : PJ_lonlat
	{
		public override string Name { get { return "longlat"; } }
		public override string DescriptionName { get { return "Lat/long (Geodetic alias)"; } }
	}
}
