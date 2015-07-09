//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Public (application) include file for PROJ.4 API, and constants.
// Author:	Frank Warmerdam, <warmerdam@pobox.com>
//
//*****************************************************************************
// Copyright (c) 2001, Frank Warmerdam <warmerdam@pobox.com>
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

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// Try to update this every version!
		public const int PJ_VERSION=491;

		public const double RAD_TO_DEG=57.29577951308232;
		public const double DEG_TO_RAD=0.0174532925199432958;

		const double EPS7=1.0e-7;
		const double EPS11=1.0e-11;
		const double EPS12=1.0e-12;
		const double TOL10=1.0e-10;

		public static int pj_errno; // global error return code // TODO wegmachen
	}

	public enum PJ_LOG : int
	{
		NONE=0,
		ERROR=1,
		DEBUG_MAJOR=2,
		DEBUG_MINOR=3,
	}
}
