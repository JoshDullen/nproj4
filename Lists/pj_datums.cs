//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Built in datum list.
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

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// The ellipse code must match one from pj_ellps.cs. The datum id should
		// be kept to 12 characters or less if possible. Use the official OGC
		// datum name for the comments if available.
		public static string GetDatum(string name, out string defn, out string ellipse_id)
		{
			switch(name)
			{
				case "WGS84": defn="towgs84=0,0,0"; ellipse_id="WGS84"; return "";
				case "GGRS87": defn="towgs84=-199.87,74.79,246.62"; ellipse_id="GRS80"; return "Greek_Geodetic_Reference_System_1987";
				case "NAD83": defn="towgs84=0,0,0"; ellipse_id="GRS80"; return "North_American_Datum_1983";
				case "NAD27": defn="nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat"; ellipse_id="clrk66"; return "North_American_Datum_1927";
				case "potsdam": defn="towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7"; ellipse_id="bessel"; return "Potsdam Rauenberg 1950 DHDN";
				case "carthage": defn="towgs84=-263.0,6.0,431.0"; ellipse_id="clrk80ign"; return "Carthage 1934 Tunisia";
				case "hermannskogel": defn="towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232"; ellipse_id="bessel"; return "Hermannskogel";
				case "ire65": defn="towgs84=482.530,-130.596,564.557,-1.042,-0.214,-0.631,8.15"; ellipse_id="mod_airy"; return "Ireland 1965";
				case "nzgd49": defn="towgs84=59.47,-5.04,187.44,0.47,-0.1,1.024,-4.5993"; ellipse_id="intl"; return "New Zealand Geodetic Datum 1949";
				case "OSGB36": defn="towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894"; ellipse_id="airy"; return "Airy 1830";
			}

			defn=ellipse_id=null;
			return null;
		}

		public static string GetPrimeMeridian(string name)
		{
			switch(name)
			{
				case "greenwich": return "0dE";
				case "lisbon": return "9d07'54.862\"W";
				case "paris": return "2d20'14.025\"E";
				case "bogota": return "74d04'51.3\"W";
				case "madrid": return "3d41'16.58\"W";
				case "rome": return "12d27'8.4\"E";
				case "bern": return "7d26'22.5\"E";
				case "jakarta": return "106d48'27.79\"E";
				case "ferro": return "17d40'W";
				case "brussels": return "4d22'4.71\"E";
				case "stockholm": return "18d3'29.8\"E";
				case "athens": return "23d42'58.815\"E";
				case "oslo": return "10d43'22.5\"E";
				default: return null;
			}
		}
	}
}
