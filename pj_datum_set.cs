//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Apply datum definition to PJ structure from initialization string.
// Author:	Frank Warmerdam, warmerda@home.com
//
//*****************************************************************************
// Copyright (c) 2000, Frank Warmerdam
// Copyright (c) 2008-2015 by the Authors
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

using System;
using System.Collections.Generic;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// SEC_TO_RAD = Pi/180/3600
		const double SEC_TO_RAD=4.84813681109535993589914102357e-6;

		//**********************************************************************
		//							pj_datum_set()
		//**********************************************************************
		static bool pj_datum_set(projCtx ctx, List<string> pl, PJ projdef)
		{
			projdef.datum_type=PJD.UNKNOWN;

			// --------------------------------------------------------------------
			//		Is there a datum definition in the parameters list? If so,
			//		add the defining values to the parameter list. Note that
			//		this will append the ellipse definition as well as the
			//		towgs84= and related parameters. It should also be pointed
			//		out that the addition is permanent rather than temporary
			//		like most other keyword expansion so that the ellipse
			//		definition will last into the pj_ell_set() function called
			//		after this one.
			// --------------------------------------------------------------------
			string name=pj_param_s(ctx, pl, "datum");
			if(name!=null&&name!="")
			{
				string defn, ellipse_id;
				if(GetDatum(name, out defn, out ellipse_id)==null) { pj_ctx_set_errno(ctx, -9); return true; }

				if(ellipse_id!=null&&ellipse_id.Length>0) pl.Add(pj_mkparam("ellps="+ellipse_id));
				if(defn!=null&&defn.Length>0) pl.Add(pj_mkparam(defn));
			}

			string nadgrids=pj_param_s(ctx, pl, "nadgrids");
			string towgs84=pj_param_s(ctx, pl, "towgs84");
			string catalog=pj_param_s(ctx, pl, "catalog");

			// --------------------------------------------------------------------
			//		Check for nadgrids parameter.
			// --------------------------------------------------------------------
			if(!string.IsNullOrEmpty(nadgrids))
			{
				// We don't actually save the value separately. It will continue
				// to exist in the param list for use in pj_apply_gridshift.cs
				projdef.datum_type=PJD.GRIDSHIFT;
			}

			// --------------------------------------------------------------------
			//		Check for grid catalog parameter, and optional date.
			// --------------------------------------------------------------------
			else if(!string.IsNullOrEmpty(catalog))
			{
				projdef.datum_type=PJD.GRIDSHIFT;
				projdef.catalog_name=catalog;

				string date=pj_param_s(ctx, pl, "date");
				if(!string.IsNullOrEmpty(date)) projdef.datum_date=Free.Ports.Proj4.Gridshift.Grid.pj_gc_parsedate(ctx, date);
			}

			// --------------------------------------------------------------------
			//		Check for towgs84 parameter.
			// --------------------------------------------------------------------
			else if(!string.IsNullOrEmpty(towgs84))
			{
				projdef.datum_params[0]=projdef.datum_params[1]=projdef.datum_params[2]=0;
				projdef.datum_params[3]=projdef.datum_params[4]=projdef.datum_params[5]=projdef.datum_params[6]=0;

				// parse out the parameters
				string s=towgs84;
				string[] ss=s.Split(',');
				for(int i=0; i<7&&i<ss.Length; i++) projdef.datum_params[i]=double.Parse(ss[i], nc);

				if(projdef.datum_params[3]!=0.0||projdef.datum_params[4]!=0.0||projdef.datum_params[5]!=0.0||projdef.datum_params[6]!=0.0)
				{
					projdef.datum_type=PJD._7PARAM;

					// transform from arc seconds to radians
					projdef.datum_params[3]*=SEC_TO_RAD;
					projdef.datum_params[4]*=SEC_TO_RAD;
					projdef.datum_params[5]*=SEC_TO_RAD;

					// transform from parts per million to scaling factor
					projdef.datum_params[6]=(projdef.datum_params[6]/1000000.0)+1;
				}
				else projdef.datum_type=PJD._3PARAM;

				// Note that pj_init() will later switch datum_type to PJD.WGS84 if shifts are all zero, and ellipsoid is WGS84 or GRS80
			}

			return false;
		}
	}
}
