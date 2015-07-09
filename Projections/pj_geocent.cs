//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Stub projection for geocentric. The transformation isn't
//			really done here since this code is 2D. The real transformation
//			is handled by pj_transform.cs.
// Author:	Frank Warmerdam, warmerdam@pobox.com
//
//*****************************************************************************
// Copyright (c) 2002, Frank Warmerdam
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

namespace Free.Ports.Proj4.Projections
{
	class PJ_geocent : PJ
	{
		public override string Name { get { return "geocent"; } }
		public override string DescriptionName { get { return "Geocentric"; } }
		public override string DescriptionType { get { return "Geoc"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		XY forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=lp.lam;
			xy.y=lp.phi;

			return xy;
		}

		LP inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=xy.y;
			lp.lam=xy.x;

			return lp;
		}

		public override PJ Init()
		{
			is_geocent=true;
			x0=0.0;
			y0=0.0;
			inv=inverse;
			fwd=forward;

			return this;
		}
	}
}
