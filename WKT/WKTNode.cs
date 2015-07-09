using System;
using System.Collections.Generic;
using System.Globalization;

namespace Free.Ports.Proj4.WKT
{
	public class WKTNode
	{
		static CultureInfo nc=new CultureInfo("");
		public string Key;
		public List<string> Values=new List<string>();
		public List<WKTNode> Children=new List<WKTNode>();

		public List<string> GetValuesOfFirstNode(string key)
		{
			if(Key==key) return Values;

			foreach(WKTNode child in Children)
			{
				List<string> ret=child.GetValuesOfFirstNode(key);
				if(ret!=null) return ret;
			}

			return null;
		}

		string GetFirstValueOfFirstNode(string key)
		{
			List<string> tmp=GetValuesOfFirstNode(key);
			if(tmp==null) return null;
			return tmp[0];
		}

		public List<WKTNode> GetNodes(string key)
		{
			List<WKTNode> ret=new List<WKTNode>();
			if(Key==key) ret.Add(this);
			foreach(WKTNode child in Children) ret.AddRange(child.GetNodes(key));
			return ret;
		}

		public bool ContainsKey(string key)
		{
			if(Key==key) return true;

			foreach(WKTNode child in Children)
			{
				bool ret=child.ContainsKey(key);
				if(ret) return true;
			}

			return false;
		}

		#region ToProj4
		const string SRS_PT_ALBERS_CONIC_EQUAL_AREA="Albers_Conic_Equal_Area";
		const string SRS_PT_AZIMUTHAL_EQUIDISTANT="Azimuthal_Equidistant";
		const string SRS_PT_BONNE="Bonne";
		const string SRS_PT_CASSINI_SOLDNER="Cassini_Soldner";
		const string SRS_PT_CYLINDRICAL_EQUAL_AREA="Cylindrical_Equal_Area";
		const string SRS_PT_ECKERT_IV="Eckert_IV";
		const string SRS_PT_ECKERT_VI="Eckert_VI";
		const string SRS_PT_EQUIDISTANT_CONIC="Equidistant_Conic";
		const string SRS_PT_EQUIRECTANGULAR="Equirectangular";
		const string SRS_PT_GALL_STEREOGRAPHIC="Gall_Stereographic";
		const string SRS_PT_GNOMONIC="Gnomonic";
		const string SRS_PT_HOTINE_OBLIQUE_MERCATOR="Hotine_Oblique_Mercator";
		const string SRS_PT_LABORDE_OBLIQUE_MERCATOR="Laborde_Oblique_Mercator";
		const string SRS_PT_LAMBERT_CONFORMAL_CONIC_1SP="Lambert_Conformal_Conic_1SP";
		const string SRS_PT_LAMBERT_CONFORMAL_CONIC_2SP="Lambert_Conformal_Conic_2SP";
		const string SRS_PT_LAMBERT_CONFORMAL_CONIC_2SP_BELGIUM="Lambert_Conformal_Conic_2SP_Belgium)";
		const string SRS_PT_LAMBERT_AZIMUTHAL_EQUAL_AREA="Lambert_Azimuthal_Equal_Area";
		const string SRS_PT_MERCATOR_1SP="Mercator_1SP";
		const string SRS_PT_MERCATOR_2SP="Mercator_2SP";
		const string SRS_PT_MILLER_CYLINDRICAL="Miller_Cylindrical";
		const string SRS_PT_MOLLWEIDE="Mollweide";
		const string SRS_PT_NEW_ZEALAND_MAP_GRID="New_Zealand_Map_Grid";
		const string SRS_PT_OBLIQUE_STEREOGRAPHIC="Oblique_Stereographic";
		const string SRS_PT_ORTHOGRAPHIC="Orthographic";
		const string SRS_PT_POLAR_STEREOGRAPHIC="Polar_Stereographic";
		const string SRS_PT_POLYCONIC="Polyconic";
		const string SRS_PT_ROBINSON="Robinson";
		const string SRS_PT_SINUSOIDAL="Sinusoidal";
		const string SRS_PT_STEREOGRAPHIC="Stereographic";
		const string SRS_PT_SWISS_OBLIQUE_CYLINDRICAL="Swiss_Oblique_Cylindrical";
		const string SRS_PT_TRANSVERSE_MERCATOR="Transverse_Mercator";
		const string SRS_PT_TRANSVERSE_MERCATOR_SOUTH_ORIENTED="Transverse_Mercator_South_Orientated";
		const string SRS_PT_TUNISIA_MINING_GRID="Tunisia_Mining_Grid";
		const string SRS_PT_VANDERGRINTEN="VanDerGrinten";

		// special mapinfo variants on Transverse Mercator
		const string SRS_PT_TRANSVERSE_MERCATOR_MI_21="Transverse_Mercator_MapInfo_21";
		const string SRS_PT_TRANSVERSE_MERCATOR_MI_22="Transverse_Mercator_MapInfo_22";
		const string SRS_PT_TRANSVERSE_MERCATOR_MI_23="Transverse_Mercator_MapInfo_23";
		const string SRS_PT_TRANSVERSE_MERCATOR_MI_24="Transverse_Mercator_MapInfo_24";
		const string SRS_PT_TRANSVERSE_MERCATOR_MI_25="Transverse_Mercator_MapInfo_25";

		const string SRS_PT_TWO_POINT_EQUIDISTANT="Two_Point_Equidistant";
		const string SRS_PT_KROVAK="Krovak";
		const string SRS_PT_HOTINE_OBLIQUE_MERCATOR_TWO_POINT_NATURAL_ORIGIN="Hotine_Oblique_Mercator_Two_Point_Natural_Origin";
		const string SRS_PT_GOODE_HOMOLOSINE="Goode_Homolosine";
		const string SRS_PT_GEOSTATIONARY_SATELLITE="Geostationary_Satellite";

		const string SRS_PP_CENTRAL_MERIDIAN="central_meridian";
		const string SRS_PP_SCALE_FACTOR="scale_factor";
		const string SRS_PP_STANDARD_PARALLEL_1="standard_parallel_1";
		const string SRS_PP_STANDARD_PARALLEL_2="standard_parallel_2";
		const string SRS_PP_LONGITUDE_OF_CENTER="longitude_of_center";
		const string SRS_PP_LATITUDE_OF_CENTER="latitude_of_center";
		const string SRS_PP_LONGITUDE_OF_ORIGIN="longitude_of_origin";
		const string SRS_PP_LATITUDE_OF_ORIGIN="latitude_of_origin";
		const string SRS_PP_FALSE_EASTING="false_easting";
		const string SRS_PP_FALSE_NORTHING="false_northing";
		const string SRS_PP_AZIMUTH="azimuth";
		const string SRS_PP_LONGITUDE_OF_POINT_1="longitude_of_point_1";
		const string SRS_PP_LATITUDE_OF_POINT_1="latitude_of_point_1";
		const string SRS_PP_LONGITUDE_OF_POINT_2="longitude_of_point_2";
		const string SRS_PP_LATITUDE_OF_POINT_2="latitude_of_point_2";
		const string SRS_PP_LONGITUDE_OF_POINT_3="longitude_of_point_3";
		const string SRS_PP_LATITUDE_OF_POINT_3="latitude_of_point_3";
		const string SRS_PP_RECTIFIED_GRID_ANGLE="rectified_grid_angle";

		// for Geostationary_Satellite
		const string SRS_PP_SATELLITE_HEIGHT="satellite_height";

		// for Two_Point_Equidistant
		const string SRS_PP_LATITUDE_OF_1ST_POINT="Latitude_Of_1st_Point";
		const string SRS_PP_LONGITUDE_OF_1ST_POINT="Longitude_Of_1st_Point";
		const string SRS_PP_LATITUDE_OF_2ND_POINT="Latitude_Of_2nd_Point";
		const string SRS_PP_LONGITUDE_OF_2ND_POINT="Longitude_Of_2nd_Point";

		const string SRS_UL_METER="Meter";
		const string SRS_UL_FOOT="Foot (International)"; // or just "FOOT"?
		const string SRS_UL_US_FOOT="U.S. Foot"; // or "US survey foot"
		const string SRS_UL_NAUTICAL_MILE="Nautical Mile";

		const string SRS_DN_NAD27="North_American_Datum_1927";
		const string SRS_DN_NAD83="North_American_Datum_1983";
		const string SRS_DN_WGS72="WGS_1972";
		const string SRS_DN_WGS84="WGS_1984";

		const double SRS_WGS84_SEMIMAJOR=6378137.0;
		const double SRS_WGS84_INVFLATTENING=298.257223563;

		const double SRS_UA_DEGREE_CONV=0.0174532925199433;

		double dfToMeter;
		double dfToDegrees;
		double dfFromGreenwich=0.0;

		/// <summary>
		/// Set the internal information for normalizing linear, and angular values.
		/// </summary>
		void GetNormInfo()
		{
			// Initialize values.
			dfFromGreenwich=GetPrimeMeridian();
			string dummy;
			dfToMeter=GetLinearUnits(out dummy);
			dfToDegrees=GetAngularUnits()/SRS_UA_DEGREE_CONV;

			if(Math.Abs(dfToDegrees-1.0)<0.000000001)
				dfToDegrees=1.0;
		}

		double GetNormProjParm(string pszName, double dfDefaultValue)
		{
			GetNormInfo();

			double dfRawResult=GetProjParm(pszName, dfDefaultValue);

			if(dfToDegrees!=1.0&&IsAngularParameter(pszName))
				return dfRawResult*dfToDegrees;

			if(dfToMeter!=1.0&&IsLinearParameter(pszName))
				return dfRawResult*dfToMeter;

			return dfRawResult;
		}

		/// <summary>
		/// Is the passed projection parameter an angular one?
		/// </summary>
		/// <param name="pszParameterName"></param>
		/// <returns></returns>
		bool IsAngularParameter(string pszParameterName)
		{
			return pszParameterName.StartsWith("long")||pszParameterName.StartsWith("lati")||
				pszParameterName==SRS_PP_CENTRAL_MERIDIAN||pszParameterName.StartsWith("standard_parallel")||
				pszParameterName==SRS_PP_AZIMUTH||pszParameterName==SRS_PP_RECTIFIED_GRID_ANGLE;
		}

		/// <summary>
		/// Is the passed projection parameter an linear one measured in
		/// meters or some similar linear measure.
		/// </summary>
		/// <param name="pszParameterName"></param>
		/// <returns></returns>
		bool IsLinearParameter(string pszParameterName)
		{
			return pszParameterName.StartsWith("false")||pszParameterName==SRS_PP_SATELLITE_HEIGHT;
		}
		
		/// <summary>
		/// Fetch a projection parameter value.
		/// 
		/// NOTE: This code should be modified to translate non degree angles into
		/// degrees based on the GEOGCS unit. This has not yet been done.
		/// </summary>
		/// <param name="pszName">the name of the parameter to fetch, from the set of SRS_PP codes in ogr_srs_api.h.</param>
		/// <param name="dfDefaultValue">the value to return if this parameter doesn't exist.</param>
		/// <returns>value of parameter.</returns>
		double GetProjParm(string pszName, double dfDefaultValue)
		{
			List<WKTNode> nodes=GetNodes("PROJCS");

			if(nodes!=null)
			{
				List<WKTNode> parameters=nodes[0].GetNodes("PARAMETER");

				foreach(WKTNode param in parameters)
				{
					if(param.Values.Count<2) continue;

					if(pszName.Equals(param.Values[0], StringComparison.OrdinalIgnoreCase)) return double.Parse(param.Values[1], nc);
				}
			}

			// Try similar names, for selected parameters.
			if(pszName==SRS_PP_LATITUDE_OF_ORIGIN) return GetProjParm(SRS_PP_LATITUDE_OF_CENTER, dfDefaultValue);

			if(pszName==SRS_PP_CENTRAL_MERIDIAN)
			{
				double dfValue=GetProjParm(SRS_PP_LONGITUDE_OF_CENTER, 0.0);
				if(dfValue!=0) return dfValue;

				dfValue=GetProjParm(SRS_PP_LONGITUDE_OF_ORIGIN, 0.0);
				if(dfValue!=0) return dfValue;
			}

			// Return default value on failure.
			return dfDefaultValue;
		}

		/// <summary>
		/// Fetch angular geographic coordinate system units.
		/// 
		/// If no units are available, a value of "degree" and SRS_UA_DEGREE_CONV
		/// will be assumed. This method only checks directly under the GEOGCS node
		/// for units.
		/// </summary>
		/// <returns>the value to multiply by angular distances to transform them to radians.</returns>
		double GetAngularUnits()
		{
			List<WKTNode> nodes=GetNodes("GEOGCS");
			if(nodes==null) return SRS_UA_DEGREE_CONV;

			foreach(WKTNode node in nodes)
			{
				List<WKTNode> units=node.GetNodes("UNIT");
				foreach(WKTNode unit in units)
				{
					if(unit.Values.Count<2) continue;

					string name=unit.Values[0];
					if(name==null) continue;

					if(name.Equals("Degree", StringComparison.OrdinalIgnoreCase)||
						name.Equals("Gon", StringComparison.OrdinalIgnoreCase)||
						name.Equals("Grad", StringComparison.OrdinalIgnoreCase)||
						name.Equals("Microradian", StringComparison.OrdinalIgnoreCase)||
						name.Equals("Mil_6400", StringComparison.OrdinalIgnoreCase)||
						name.Equals("Minute", StringComparison.OrdinalIgnoreCase)||
						name.Equals("Minute_Centesimal", StringComparison.OrdinalIgnoreCase)||
						name.Equals("Radian", StringComparison.OrdinalIgnoreCase)||
						name.Equals("Second", StringComparison.OrdinalIgnoreCase)||
						name.Equals("Second_Centesimal", StringComparison.OrdinalIgnoreCase)) return double.Parse(unit.Values[1], nc);
				}
			}

			return SRS_UA_DEGREE_CONV;
		}

		/// <summary>
		/// Fetch linear projection units.
		///
		/// If no units are available, a value of "Meters" and 1.0 will be assumed.
		/// This method only checks directly under the PROJCS or LOCAL_CS node for
		/// units.
		/// </summary>
		/// <param name="ppszName">The returned value remains internal to the OGRSpatialReference
		/// and shouldn't be freed, or modified. It may be invalidated on the next
		/// OGRSpatialReference call.</param>
		/// <returns>the value to multiply by linear distances to transform them to
		/// meters.</returns>
		double GetLinearUnits(out string ppszName)
		{
			List<WKTNode> nodes=GetNodes("PROJCS");
			if(nodes==null) nodes=GetNodes("LOCAL_CS");

			ppszName="unknown";
			if(nodes==null) return 1.0;

			foreach(WKTNode node in nodes)
			{
				List<WKTNode> units=node.GetNodes("UNIT");
				foreach(WKTNode unit in units)
				{
					if(unit.Values.Count<2) continue;

					string name=unit.Values[0];
					if(name!=null||name.Length>0)
					{
						if(name.Equals("Degree", StringComparison.OrdinalIgnoreCase)||
							name.Equals("Gon", StringComparison.OrdinalIgnoreCase)||
							name.Equals("Grad", StringComparison.OrdinalIgnoreCase)||
							name.Equals("Microradian", StringComparison.OrdinalIgnoreCase)||
							name.Equals("Mil_6400", StringComparison.OrdinalIgnoreCase)||
							name.Equals("Minute", StringComparison.OrdinalIgnoreCase)||
							name.Equals("Minute_Centesimal", StringComparison.OrdinalIgnoreCase)||
							name.Equals("Radian", StringComparison.OrdinalIgnoreCase)||
							name.Equals("Second", StringComparison.OrdinalIgnoreCase)||
							name.Equals("Second_Centesimal", StringComparison.OrdinalIgnoreCase)) continue;

						ppszName=name;
					}

					return double.Parse(unit.Values[1], nc);
				}
			}

			return 1.0;
		}

		/// <summary>
		/// Fetch prime meridian info.
		/// 
		/// Returns the offset of the prime meridian from greenwich in degrees,
		/// and the prime meridian name (if requested). If no PRIMEM value exists
		/// in the coordinate system definition a value of "Greenwich" and an
		/// offset of 0.0 is assumed.
		/// 
		/// If the prime meridian name is returned, the pointer is to an internal
		/// copy of the name. It should not be freed, altered or depended on after
		/// the next OGR call.
		/// </summary>
		/// <returns>the offset to the GEOGCS prime meridian from greenwich in decimal degrees.</returns>
		double GetPrimeMeridian()
		{
			List<WKTNode> nodes=GetNodes("PRIMEM");

			if(nodes!=null)
			{
				WKTNode primem=nodes[0];
				if(primem.Values.Count>=2) return double.Parse(primem.Values[1], nc);
			}

			return 0.0;
		}

		/// <summary>
		/// Get utm zone information.
		/// 
		/// This is the same as the C function OSRGetUTMZone().
		/// </summary>
		/// <param name="pbNorth">set to true if northern hemisphere, or false if southern.</param>
		/// <returns>UTM zone number or zero if this isn't a UTM definition.</returns>
		int GetUTMZone(out bool pbNorth)
		{
			pbNorth=true;

			string pszProjection=GetFirstValueOfFirstNode("PROJECTION");

			if(pszProjection==null||pszProjection!=SRS_PT_TRANSVERSE_MERCATOR) return 0;

			if(GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0)!=0.0) return 0;
			if(GetProjParm(SRS_PP_SCALE_FACTOR, 1.0)!=0.9996) return 0;

			if(Math.Abs(GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0)-500000.0)>0.001) return 0;

			double dfFalseNorthing=GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0);
			if(dfFalseNorthing!=0.0&&Math.Abs(dfFalseNorthing-10000000.0)>0.001) return 0;

			pbNorth=dfFalseNorthing==0;

			double dfCentralMeridian=GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0);
			double dfZone=(dfCentralMeridian+183)/6.0+0.000000001;

			if(Math.Abs(dfZone-(int)dfZone)>0.00001||dfCentralMeridian<-177.00001||dfCentralMeridian>177.000001) return 0;

			return (int)dfZone;
		}

		/// <summary>
		/// Get spheroid semi major axis.
		/// </summary>
		/// <returns>semi-major axis, or SRS_WGS84_SEMIMAJOR if it can't be found.</returns>
		double GetSemiMajor()
		{
			List<string> spheroidValues=GetValuesOfFirstNode("SPHEROID");

			if(spheroidValues!=null&&spheroidValues.Count>=3)
				return double.Parse(spheroidValues[1], nc);

			return SRS_WGS84_SEMIMAJOR;
		}

		/// <summary>
		/// Get spheroid semi minor axis.
		/// </summary>
		/// <returns>semi-minor axis, or WGS84 semi minor if it can't be found.</returns>
		double GetSemiMinor()
		{
			double dfSemiMajor=GetSemiMajor();
			double dfInvFlattening=GetInvFlattening();

			if(Math.Abs(dfInvFlattening)<0.000000000001) return dfSemiMajor;
			return dfSemiMajor*(1.0-1.0/dfInvFlattening);
		}

		/// <summary>
		/// Get spheroid inverse flattening.
		/// </summary>
		/// <returns>inverse flattening, or SRS_WGS84_INVFLATTENING if it can't be found.</returns>
		double GetInvFlattening()
		{
			List<string> spheroidValues=GetValuesOfFirstNode("SPHEROID");

			if(spheroidValues!=null&&spheroidValues.Count>=3)
				return double.Parse(spheroidValues[2], nc);

			return SRS_WGS84_INVFLATTENING;
		}

		/// <summary>
		/// Get the authority name for a node.
		///
		/// This method is used to query an AUTHORITY[] node from within the 
		/// WKT tree, and fetch the authority name value.
		///
		/// The most common authority is "EPSG".
		/// </summary>
		/// <param name="pszTargetKey"> the partial or complete path to the node to get an authority from. ie. "PROJCS", "GEOGCS", "GEOGCS|UNIT" or NULL to search for an authority node on the root element.</param>
		/// <returns>value code from authority node, or NULL on failure. The value returned is internal and should not be freed or modified.</returns>
		string GetAuthorityName(string pszTargetKey)
		{
			List<WKTNode> nodes=GetNodes(pszTargetKey);

			if(nodes!=null)
			{
				WKTNode target=nodes[0];

				List<WKTNode> authorities=target.GetNodes("AUTHORITY");

				if(authorities.Count>0)
				{
					WKTNode autority=authorities[0];
					return autority.Values[0];
				}
			}

			return null;
		}

		/// <summary>
		/// Get the authority code for a node.
		///
		/// This method is used to query an AUTHORITY[] node from within the 
		/// WKT tree, and fetch the code value.  
		///
		/// While in theory values may be non-numeric, for the EPSG authority all
		/// code values should be integral.
		/// </summary>
		/// <param name="pszTargetKey">the partial or complete path to the node to
		/// get an authority from.  ie. "PROJCS", "GEOGCS", "GEOGCS|UNIT" or NULL to
		/// search for an authority node on the root element.</param>
		/// <returns>value code from authority node, or NULL on failure. The value
		/// returned is internal and should not be freed or modified.</returns>
		string GetAuthorityCode(string pszTargetKey)
		{
			List<WKTNode> nodes=GetNodes(pszTargetKey);

			if(nodes!=null)
			{
				WKTNode target=nodes[0];

				List<WKTNode> authorities=target.GetNodes("AUTHORITY");

				if(authorities!=null)
				{
					WKTNode autority=authorities[0];
					return autority.Values[1];
				}
			}

			return null;
		}

		// ***********************************************************************
		//						EPSGGetWGS84Transform()
		//
		//		The following code attempts to find a bursa-wolf
		//		transformation from this GeogCS to WGS84 (4326).
		//
		//		Faults:
		//		*	I think there are codes other than 9603 and 9607 that
		//			return compatible, or easily transformed parameters.
		//		*	Only the first path from the given GeogCS is checked due
		//			to limitations in the CSV API.
		// ***********************************************************************
		bool EPSGGetWGS84Transform(int nGeogCS, out double[] padfTransform)
		{
			// Out Variable definieren
			padfTransform=new double[7];

			string szCode=string.Format("{0}", nGeogCS);

			if(!coordCSV.ContainsKey(nGeogCS)) return false;

			List<string> papszLine=coordCSV[nGeogCS];

			// Verify that the method code is one of our accepted ones.
			int nMethodCode=int.Parse(papszLine[10]); //"COORD_OP_METHOD_CODE"

			if(nMethodCode!=9603&&nMethodCode!=9607&&nMethodCode!=9606)
				return false;

			// Fetch the transformation parameters.
			int iDXField=11; //"DX"

			for(int iField=0; iField<7; iField++)
			{
				padfTransform[iField]=double.Parse(papszLine[iDXField+iField], nc);
			}

			// --------------------------------------------------------------------
			// 9607 - coordinate frame rotation has reverse signs on the
			// rotational coefficients.  Fix up now since we internal
			// operate according to method 9606 (position vector 7-parameter).
			// --------------------------------------------------------------------
			if(nMethodCode==9607)
			{
				padfTransform[3]*=-1;
				padfTransform[4]*=-1;
				padfTransform[5]*=-1;
			}

			return true;
		}

		public string ToProj4String()
		{
			// Werte ermitteln
			string PROJECTION=GetFirstValueOfFirstNode("PROJECTION");
			List<string> PRIMEM=GetValuesOfFirstNode("PRIMEM");
			bool geographic=ContainsKey("GEOGCS");

			#region Handle the projection definition.
			if(PROJECTION==null&&!geographic) return ""; // LOCAL_CS, or incompletely initialized coordinate systems.

			string szProj4="";
			if(PROJECTION==null&&geographic)
			{
				szProj4="+proj=longlat ";
			}
			else if(PROJECTION==SRS_PT_CYLINDRICAL_EQUAL_AREA)
			{
				szProj4=string.Format(nc, "+proj=cea +lon_0={0:F16} +lat_ts{1:F16} +x_0={2:F16} +y_0={3:F16} ",
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0), GetNormProjParm(SRS_PP_STANDARD_PARALLEL_1, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0), GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_BONNE)
			{
				szProj4=string.Format(nc, "+proj=bonne +lon_0={0:F16} +lat_1={1:F16} +x_0={2:F16} +y_0={3:F16} ",
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_STANDARD_PARALLEL_1, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_CASSINI_SOLDNER)
			{
				szProj4=string.Format(nc, "+proj=cass +lat_0={0:F16} +lon_0={1:F16} +x_0={2:F16} +y_0={3:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_NEW_ZEALAND_MAP_GRID)
			{
				szProj4=string.Format(nc, "+proj=nzmg +lat_0={0:F16} +lon_0={1:F16} +x_0={2:F16} +y_0={3:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_TRANSVERSE_MERCATOR||
					PROJECTION==SRS_PT_TRANSVERSE_MERCATOR_MI_21||
					PROJECTION==SRS_PT_TRANSVERSE_MERCATOR_MI_22||
					PROJECTION==SRS_PT_TRANSVERSE_MERCATOR_MI_23||
					PROJECTION==SRS_PT_TRANSVERSE_MERCATOR_MI_24||
					PROJECTION==SRS_PT_TRANSVERSE_MERCATOR_MI_25)
			{
				bool bNorth=false;
				int nZone=GetUTMZone(out bNorth);

				if(nZone!=0)
				{
					szProj4=string.Format(nc, "+proj=utm +zone={0} ", nZone);
					if(!bNorth) szProj4+="+south ";
				}
				else
				{
					szProj4=string.Format(nc, "+proj=tmerc +lat_0={0:F16} +lon_0={1:F16} +k={2:F16} +x_0={3:F16} +y_0={4:F16} ",
						GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
						GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
						GetNormProjParm(SRS_PP_SCALE_FACTOR, 1.0),
						GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
						GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
				}
			}
			else if(PROJECTION==SRS_PT_MERCATOR_1SP)
			{
				szProj4=string.Format(nc, "+proj=merc +lon_0={0:F16} +k={1:F16} +x_0={2:F16} +y_0={3:F16} ",
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_SCALE_FACTOR, 1.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_MERCATOR_2SP)
			{
				szProj4=string.Format(nc, "+proj=merc +lon_0={0:F16} +lat_ts={1:F16} +x_0={2:F16} +y_0={3:F16} ",
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_STANDARD_PARALLEL_1, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_OBLIQUE_STEREOGRAPHIC)
			{
				szProj4=string.Format(nc, "+proj=sterea +lat_0={0:F16} +lon_0={1:F16} +k={2:F16} +x_0={3:F16} +y_0={4:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_SCALE_FACTOR, 1.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_STEREOGRAPHIC)
			{
				szProj4=string.Format(nc, "+proj=stere +lat_0={0:F16} +lon_0={1:F16} +k={2:F16} +x_0={3:F16} +y_0={4:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_SCALE_FACTOR, 1.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_POLAR_STEREOGRAPHIC)
			{
				if(GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0)>=0.0)
				{
					szProj4=string.Format(nc, "+proj=stere +lat_0=90 +lat_ts={0:F16} +lon_0={1:F16} +k={2:F16} +x_0={3:F16} +y_0={4:F16} ",
						GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 90.0),
						GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
						GetNormProjParm(SRS_PP_SCALE_FACTOR, 1.0),
						GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
						GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
				}
				else
				{
					szProj4=string.Format(nc, "+proj=stere +lat_0=-90 +lat_ts={0:F16} +lon_0={1:F16} +k={2:F16} +x_0={3:F16} +y_0={4:F16} ",
						GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, -90.0),
						GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
						GetNormProjParm(SRS_PP_SCALE_FACTOR, 1.0),
						GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
						GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
				}
			}
			else if(PROJECTION==SRS_PT_EQUIRECTANGULAR)
			{
				szProj4=string.Format(nc, "+proj=eqc +lat_ts={0:F16} +lon_0={1:F16} +x_0={2:F16} +y_0={3:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_GNOMONIC)
			{
				szProj4=string.Format(nc, "+proj=gnom +lat_0={0:F16} +lon_0={1:F16} +x_0={2:F16} +y_0={3:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_ORTHOGRAPHIC)
			{
				szProj4=string.Format(nc, "+proj=ortho +lat_0={0:F16} +lon_0={1:F16} +x_0={2:F16} +y_0={3:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_LAMBERT_AZIMUTHAL_EQUAL_AREA)
			{
				szProj4=string.Format(nc, "+proj=laea +lat_0={0:F16} +lon_0={1:F16} +x_0={2:F16} +y_0={3:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_AZIMUTHAL_EQUIDISTANT)
			{
				szProj4=string.Format(nc, "+proj=aeqd +lat_0={0:F16} +lon_0={1:F16} +x_0={2:F16} +y_0={3:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_EQUIDISTANT_CONIC)
			{
				szProj4=string.Format(nc, "+proj=eqdc +lat_0={0:F16} +lon_0={1:F16} +lat_1={2:F16} +lat_2={3:F16} +x_0={4:F16} +y_0={5:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_CENTER, 0.0),
					GetNormProjParm(SRS_PP_LONGITUDE_OF_CENTER, 0.0),
					GetNormProjParm(SRS_PP_STANDARD_PARALLEL_1, 0.0),
					GetNormProjParm(SRS_PP_STANDARD_PARALLEL_2, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_MILLER_CYLINDRICAL)
			{
				szProj4=string.Format(nc, "+proj=mill +lat_0={0:F16} +lon_0={1:F16} +x_0={2:F16} +y_0={3:F16} +R_A ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_MOLLWEIDE)
			{
				szProj4=string.Format(nc, "+proj=moll +lon_0={0:F16} +x_0={1:F16} +y_0={2:F16} ",
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_ECKERT_IV)
			{
				szProj4=string.Format(nc, "+proj=eck4 +lon_0={0:F16} +x_0={1:F16} +y_0={2:F16} ",
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_ECKERT_VI)
			{
				szProj4=string.Format(nc, "+proj=eck6 +lon_0={0:F16} +x_0={1:F16} +y_0={2:F16} ",
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_POLYCONIC)
			{
				szProj4=string.Format(nc, "+proj=poly +lat_0={0:F16} +lon_0={1:F16} +x_0={2:F16} +y_0={3:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_ALBERS_CONIC_EQUAL_AREA)
			{
				szProj4=string.Format(nc, "+proj=aea +lat_1={0:F16} +lat_2={1:F16} +lat_0={2:F16} +lon_0={3:F16} +x_0={4:F16} +y_0={5:F16} ",
					GetNormProjParm(SRS_PP_STANDARD_PARALLEL_1, 0.0),
					GetNormProjParm(SRS_PP_STANDARD_PARALLEL_2, 0.0),
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_ROBINSON)
			{
				szProj4=string.Format(nc, "+proj=robin +lon_0={0:F16} +x_0={1:F16} +y_0={2:F16} ",
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_VANDERGRINTEN)
			{
				szProj4=string.Format(nc, "+proj=vandg +lon_0={0:F16} +x_0={1:F16} +y_0={2:F16} +R_A ",
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_SINUSOIDAL)
			{
				szProj4=string.Format(nc, "+proj=sinu +lon_0={0:F16} +x_0={1:F16} +y_0={2:F16} ",
					GetNormProjParm(SRS_PP_LONGITUDE_OF_CENTER, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_GALL_STEREOGRAPHIC)
			{
				szProj4=string.Format(nc, "+proj=gall +lon_0={0:F16} +x_0={1:F16} +y_0={2:F16} ",
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_GOODE_HOMOLOSINE)
			{
				szProj4=string.Format(nc, "+proj=goode +lon_0={0:F16} +x_0={1:F16} +y_0={2:F16} ",
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_GEOSTATIONARY_SATELLITE)
			{
				szProj4=string.Format(nc, "+proj=geos +lon_0={0:F16} +h={1:F16} +x_0={2:F16} +y_0={3:F16} ",
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_SATELLITE_HEIGHT, 35785831.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_LAMBERT_CONFORMAL_CONIC_2SP||PROJECTION==SRS_PT_LAMBERT_CONFORMAL_CONIC_2SP_BELGIUM)
			{
				szProj4=string.Format(nc, "+proj=lcc +lat_1={0:F16} +lat_2={1:F16} +lat_0={2:F16} +lon_0={3:F16} +x_0={4:F16} +y_0={5:F16} ",
					GetNormProjParm(SRS_PP_STANDARD_PARALLEL_1, 0.0),
					GetNormProjParm(SRS_PP_STANDARD_PARALLEL_2, 0.0),
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_LAMBERT_CONFORMAL_CONIC_1SP)
			{
				szProj4=string.Format(nc, "+proj=lcc +lat_1={0:F16} +lat_0={1:F16} +lon_0={2:F16} +k_0={3:F16} +x_0={4:F16} +y_0={5:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_SCALE_FACTOR, 1.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_HOTINE_OBLIQUE_MERCATOR)
			{
				// not clear how ProjParm[3] - angle from rectified to skewed grid -
				// should be applied ... see the +not_rot flag for PROJ.4.
				// Just ignoring for now.

				// special case for swiss oblique mercator : see bug 423
				if(Math.Abs(GetNormProjParm(SRS_PP_AZIMUTH, 0.0)-90.0)<0.0001&&Math.Abs(GetNormProjParm(SRS_PP_RECTIFIED_GRID_ANGLE, 0.0)-90.0)<0.0001)
				{
					szProj4=string.Format(nc, "+proj=somerc +lat_0={0:F16} +lon_0={1:F16} +x_0={2:F16} +y_0={3:F16} ",
						GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
						GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
						GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
						GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
				}
				else
				{
					szProj4=string.Format(nc, "+proj=omerc +lat_0={0:F16} +lonc={1:F16} +alpha={2:F16} +k={3:F16} +x_0={4:F16} +y_0={5:F16} ",
						GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
						GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
						GetNormProjParm(SRS_PP_AZIMUTH, 0.0),
						GetNormProjParm(SRS_PP_SCALE_FACTOR, 1.0),
						GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
						GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
				}
			}
			else if(PROJECTION==SRS_PT_HOTINE_OBLIQUE_MERCATOR_TWO_POINT_NATURAL_ORIGIN)
			{
				szProj4=string.Format(nc, "+proj=omerc +lat_0={0:F16} +lon_1={1:F16} +lat_1={2:F16} +lon_2={3:F16} +lat_2={4:F16} +k={5:F16} +x_0={6:F16} +y_0={7:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_LATITUDE_OF_POINT_1, 0.0),
					GetNormProjParm(SRS_PP_LONGITUDE_OF_POINT_1, 0.0),
					GetNormProjParm(SRS_PP_LATITUDE_OF_POINT_2, 0.0),
					GetNormProjParm(SRS_PP_LONGITUDE_OF_POINT_2, 0.0),
					GetNormProjParm(SRS_PP_SCALE_FACTOR, 1.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_KROVAK)
			{
				szProj4=string.Format(nc, "+proj=krovak +lat_0={0:F16} +lon_0={1:F16} +alpha={2:F16} +k={3:F16} +x_0={4:F16} +y_0={5:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_CENTER, 0.0),
					GetNormProjParm(SRS_PP_LONGITUDE_OF_CENTER, 0.0),
					GetNormProjParm(SRS_PP_AZIMUTH, 0.0),
					GetNormProjParm(SRS_PP_SCALE_FACTOR, 1.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_TWO_POINT_EQUIDISTANT)
			{
				szProj4=string.Format(nc, "+proj=tpeqd +lat_1={0:F16} +lon_1={1:F16} +lat_2={2:F16} +lon_2={3:F16} +x_0={4:F16} +y_0={5:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_1ST_POINT, 0.0),
					GetNormProjParm(SRS_PP_LONGITUDE_OF_1ST_POINT, 0.0),
					GetNormProjParm(SRS_PP_LATITUDE_OF_2ND_POINT, 0.0),
					GetNormProjParm(SRS_PP_LONGITUDE_OF_2ND_POINT, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else if(PROJECTION==SRS_PT_SWISS_OBLIQUE_CYLINDRICAL) // Note: This never really gets used currently. See bug 423
			{
				szProj4=string.Format(nc, "+proj=somerc +lat_0={0:F16} +lon_0={1:F16} +x_0={2:F16} +y_0={3:F16} ",
					GetNormProjParm(SRS_PP_LATITUDE_OF_ORIGIN, 0.0),
					GetNormProjParm(SRS_PP_CENTRAL_MERIDIAN, 0.0),
					GetNormProjParm(SRS_PP_FALSE_EASTING, 0.0),
					GetNormProjParm(SRS_PP_FALSE_NORTHING, 0.0));
			}
			else return ""; // Not supported coordinate system
			#endregion

			#region Handle earth model.  For now we just always emit the user defined ellipsoid parameters.
			double dfSemiMajor=GetSemiMajor();
			double dfInvFlattening=GetInvFlattening();
			string pszPROJ4Ellipse=null;
			string DATUM=GetFirstValueOfFirstNode("DATUM");

			if(Math.Abs(dfSemiMajor-6378249.145)<0.01&&Math.Abs(dfInvFlattening-293.465)<0.0001) pszPROJ4Ellipse="clrk80"; // Clark 1880
			else if(Math.Abs(dfSemiMajor-6378245.0)<0.01&&Math.Abs(dfInvFlattening-298.3)<0.0001) pszPROJ4Ellipse="krass"; // Krassovsky
			else if(Math.Abs(dfSemiMajor-6378388.0)<0.01&&Math.Abs(dfInvFlattening-297.0)<0.0001) pszPROJ4Ellipse="intl"; // International 1924
			else if(Math.Abs(dfSemiMajor-6378160.0)<0.01&&Math.Abs(dfInvFlattening-298.25)<0.0001) pszPROJ4Ellipse="aust_SA"; // Australian
			else if(Math.Abs(dfSemiMajor-6377397.155)<0.01&&Math.Abs(dfInvFlattening-299.1528128)<0.0001) pszPROJ4Ellipse="bessel"; // Bessel 1841
			else if(Math.Abs(dfSemiMajor-6377483.865)<0.01&&Math.Abs(dfInvFlattening-299.1528128)<0.0001) pszPROJ4Ellipse="bess_nam"; // Bessel 1841 (Namibia / Schwarzeck)
			else if(Math.Abs(dfSemiMajor-6378160.0)<0.01&&Math.Abs(dfInvFlattening-298.247167427)<0.0001) pszPROJ4Ellipse="GRS67"; // GRS 1967
			else if(Math.Abs(dfSemiMajor-6378137)<0.01&&Math.Abs(dfInvFlattening-298.257222101)<0.000001) pszPROJ4Ellipse="GRS80"; // GRS 1980
			else if(Math.Abs(dfSemiMajor-6378206.4)<0.01&&Math.Abs(dfInvFlattening-294.9786982)<0.0001) pszPROJ4Ellipse="clrk66"; // Clarke 1866
			else if(Math.Abs(dfSemiMajor-6378206.4)<0.01&&Math.Abs(dfInvFlattening-294.9786982)<0.0001) pszPROJ4Ellipse="mod_airy"; // Modified Airy
			else if(Math.Abs(dfSemiMajor-6377563.396)<0.01&&Math.Abs(dfInvFlattening-299.3249646)<0.0001) pszPROJ4Ellipse="airy"; // Modified Airy
			else if(Math.Abs(dfSemiMajor-6378200)<0.01&&Math.Abs(dfInvFlattening-298.3)<0.0001) pszPROJ4Ellipse="helmert"; // Helmert 1906
			else if(Math.Abs(dfSemiMajor-6378155)<0.01&&Math.Abs(dfInvFlattening-298.3)<0.0001) pszPROJ4Ellipse="fschr60m"; // Modified Fischer 1960
			else if(Math.Abs(dfSemiMajor-6377298.556)<0.01&&Math.Abs(dfInvFlattening-300.8017)<0.0001) pszPROJ4Ellipse="evrstSS"; // Everest (Sabah & Sarawak)
			else if(Math.Abs(dfSemiMajor-6378165.0)<0.01&&Math.Abs(dfInvFlattening-298.3)<0.0001) pszPROJ4Ellipse="WGS60";
			else if(Math.Abs(dfSemiMajor-6378145.0)<0.01&&Math.Abs(dfInvFlattening-298.25)<0.0001) pszPROJ4Ellipse="WGS66";
			else if(Math.Abs(dfSemiMajor-6378135.0)<0.01&&Math.Abs(dfInvFlattening-298.26)<0.0001) pszPROJ4Ellipse="WGS72";
			else if(Math.Abs(dfSemiMajor-6378137.0)<0.01&&Math.Abs(dfInvFlattening-298.257223563)<0.000001) pszPROJ4Ellipse="WGS84";
			else if(DATUM=="North_American_Datum_1927") pszPROJ4Ellipse="clrk66";
			else if(DATUM=="North_American_Datum_1983") pszPROJ4Ellipse="GRS80";

			if(pszPROJ4Ellipse==null) szProj4+=string.Format(nc, "+a={0:G16} +b={1:G16} ", GetSemiMajor(), GetSemiMinor());
			else szProj4+=string.Format("+ellps={0} ", pszPROJ4Ellipse);
			#endregion

			#region Translate the datum.
			string pszAuthority=GetAuthorityName("DATUM");
			int nEPSGDatum=-1;
			if(pszAuthority!=null&&pszAuthority=="EPSG")
				nEPSGDatum=int.Parse(GetAuthorityCode("DATUM"));

			string pszGeogCSAuthority=GetAuthorityName("GEOGCS");
			int nEPSGGeogCS=-1;
			if(pszGeogCSAuthority!=null&&pszGeogCSAuthority=="EPSG")
				nEPSGGeogCS=int.Parse(GetAuthorityCode("GEOGCS"));

			List<string> poTOWGS84=GetValuesOfFirstNode("TOWGS84");
			string pszPROJ4Datum=null;
			if(DATUM==null) { } // nothing
			else if(DATUM==SRS_DN_NAD27||nEPSGDatum==6267) pszPROJ4Datum="+datum=NAD27 ";
			else if(DATUM==SRS_DN_NAD83||nEPSGDatum==6269) pszPROJ4Datum="+datum=NAD83 ";
			else if(DATUM==SRS_DN_WGS84||nEPSGDatum==6326) pszPROJ4Datum="+datum=WGS84 ";
			else if(nEPSGDatum==6314) pszPROJ4Datum="+datum=potsdam ";
			else if(nEPSGDatum==6272) pszPROJ4Datum="+datum=nzgd49 ";
			else if(nEPSGDatum==6277) pszPROJ4Datum="+datum=OSGB36 ";
			else if(DATUM.Equals("wgs84", StringComparison.OrdinalIgnoreCase)||DATUM.Equals("wgs 84", StringComparison.OrdinalIgnoreCase)||
				DATUM.Equals("wgs_1984", StringComparison.OrdinalIgnoreCase)||DATUM.Equals("d_wgs_1984", StringComparison.OrdinalIgnoreCase)) pszPROJ4Datum="+datum=WGS84 ";
			else if(poTOWGS84!=null)
			{
				if(poTOWGS84.Count>2&&poTOWGS84.Count<6)
					pszPROJ4Datum=string.Format("+towgs84={0},{1},{2} ", poTOWGS84[0], poTOWGS84[1], poTOWGS84[2]);
				else if(poTOWGS84.Count>6)
					pszPROJ4Datum=string.Format("+towgs84={0},{1},{2},{3},{4},{5},{6} ", poTOWGS84[0], poTOWGS84[1], poTOWGS84[2], poTOWGS84[3], poTOWGS84[4], poTOWGS84[5], poTOWGS84[6]);
			}
			else if(nEPSGGeogCS!=-1)
			{
				double[] padfTransform;
				if(EPSGGetWGS84Transform(nEPSGGeogCS, out padfTransform))
					pszPROJ4Datum=string.Format(nc, "+towgs84={0},{1},{2},{3},{4},{5},{6} ", padfTransform[0], padfTransform[1], padfTransform[2], padfTransform[3], padfTransform[4], padfTransform[5], padfTransform[6]);
			}

			if(pszPROJ4Datum!=null)
				szProj4+=pszPROJ4Datum;
			#endregion

			#region Is there prime meridian info to apply?
			if(PRIMEM!=null&&PRIMEM.Count>=2&&double.Parse(PRIMEM[1], nc)!=0.0)
			{
				string apszAuthority=GetAuthorityName("PRIMEM");
				int nCode=-1;
				if(apszAuthority!=null&&apszAuthority=="EPSG")
				{
					nCode=int.Parse(GetAuthorityCode("PRIMEM"));
				}

				string szPMValue;
				switch(nCode)
				{
					case 8902: szPMValue="lisbon"; break;
					case 8903: szPMValue="paris"; break;
					case 8904: szPMValue="bogota"; break;
					case 8905: szPMValue="madrid"; break;
					case 8906: szPMValue="rome"; break;
					case 8907: szPMValue="bern"; break;
					case 8908: szPMValue="jakarta"; break;
					case 8909: szPMValue="ferro"; break;
					case 8910: szPMValue="brussels"; break;
					case 8911: szPMValue="stockholm"; break;
					case 8912: szPMValue="athens"; break;
					case 8913: szPMValue="oslo"; break;
					default: szPMValue=string.Format(nc, "{0:F16}", dfFromGreenwich); break;
				}

				szProj4+=string.Format("+pm={0} ", szPMValue);
			}
			#endregion

			#region Handle linear units.
			string pszLinearUnits=null;
			double dfLinearConv=GetLinearUnits(out pszLinearUnits);

			string pszPROJ4Units=null;
			if(szProj4.IndexOf("longlat")!=-1) pszPROJ4Units=null;
			else if(dfLinearConv==1.0) pszPROJ4Units="m";
			else if(dfLinearConv==1000.0) pszPROJ4Units="km";
			else if(dfLinearConv==0.0254) pszPROJ4Units="in";
			else if(pszLinearUnits==SRS_UL_FOOT) pszPROJ4Units="ft";
			else if(pszLinearUnits=="IYARD"||dfLinearConv==0.9144) pszPROJ4Units="yd";
			else if(dfLinearConv==0.001) pszPROJ4Units="mm";
			else if(dfLinearConv==0.01) pszPROJ4Units="cm";
			else if(pszLinearUnits==SRS_UL_US_FOOT) pszPROJ4Units="us-ft";
			else if(pszLinearUnits==SRS_UL_NAUTICAL_MILE) pszPROJ4Units="kmi";
			else if(pszLinearUnits=="Mile"||pszLinearUnits=="IMILE") pszPROJ4Units="mi";
			else szProj4+=string.Format(nc, "+to_meter={0:F16} ", dfLinearConv);

			if(pszPROJ4Units!=null) szProj4+=string.Format("+units={0} ", pszPROJ4Units);
			#endregion

			//Add the no_defs flag to ensure that no values from
			//proj_def.dat are implicitly used with our definitions.
			szProj4+="+no_defs ";

			return szProj4;
		}
		#endregion

		static List<string> GetCoordCSV(string line)
		{
			string[] parts=line.Split(new char[] { ',' }, StringSplitOptions.None);
			return new List<string>(parts);
		}
		
		static Dictionary<int, List<string>> coordCSV;

		static WKTNode()
		{
			coordCSV=new Dictionary<int, List<string>>();

			#region Build CSV file
			//"COORD_REF_SYS_CODE","COORD_REF_SYS_NAME","DATUM_CODE","DATUM_NAME","GREENWICH_DATUM","UOM_CODE","ELLIPSOID_CODE","PRIME_MERIDIAN_CODE","SHOW_CRS","DEPRECATED","COORD_OP_METHOD_CODE","DX","DY","DZ","RX","RY","RZ","DS"
			coordCSV.Add(4272, GetCoordCSV("4272,NZGD49,6272,New Zealand Geodetic Datum 1949,6272,9122,7022,8901,1,0,9606,59.47,-5.04,187.44,0.47,-0.1,1.024,-4.5993"));
			coordCSV.Add(4312, GetCoordCSV("4312,MGI,6312,Militar-Geographische Institut,6312,9108,7004,8901,1,0,9606,577.326,90.129,463.919,5.137,1.474,5.297,2.4232"));
			coordCSV.Add(4218, GetCoordCSV("4218,Bogota 1975,6218,Bogota 1975,6218,9122,7022,8901,1,0,9603,307,304,-318,,,,"));
			coordCSV.Add(4149, GetCoordCSV("4149,CH1903,6149,CH1903,6149,9122,7004,8901,1,0,9603,674.374,15.056,405.346,,,,"));
			coordCSV.Add(4313, GetCoordCSV("4313,Belge 1972,6313,Reseau National Belge 1972,6313,9122,7022,8901,1,0,9606,106.868628,-52.297783,103.723893,-0.336570,0.456955,-1.842183,1.0000012747"));
			coordCSV.Add(4123, GetCoordCSV("4123,KKJ,6123,Kartastokoordinaattijarjestelma 1966,6123,9122,7022,8901,1,0,9606,-96.0617,-82.4278,-121.7435,4.80107,0.34543,-1.37646,1.4964"));////////////////////////////////////////////////");

			coordCSV.Add(4001, GetCoordCSV("4001,Unknown datum based upon the Airy 1830 ellipsoid,6001,\"Not specified (based on Airy 1830 ellipsoid)\",6001,9122,7001,8901,0,0,,,,,,,,"));
			coordCSV.Add(4002, GetCoordCSV("4002,Unknown datum based upon the Airy Modified 1849 ellipsoid,6002,\"Not specified (based on Airy Modified 1849 ellipsoid)\",6002,9122,7002,8901,0,0,,,,,,,,"));
			coordCSV.Add(4003, GetCoordCSV("4003,Unknown datum based upon the Australian National Spheroid,6003,\"Not specified (based on Australian National Spheroid)\",6003,9122,7003,8901,0,0,,,,,,,,"));
			coordCSV.Add(4004, GetCoordCSV("4004,Unknown datum based upon the Bessel 1841 ellipsoid,6004,\"Not specified (based on Bessel 1841 ellipsoid)\",6004,9122,7004,8901,0,0,,,,,,,,"));
			coordCSV.Add(4005, GetCoordCSV("4005,Unknown datum based upon the Bessel Modified ellipsoid,6005,\"Not specified (based on Bessel Modified ellipsoid)\",6005,9122,7005,8901,0,0,,,,,,,,"));
			coordCSV.Add(4006, GetCoordCSV("4006,Unknown datum based upon the Bessel Namibia ellipsoid,6006,\"Not specified (based on Bessel Namibia ellipsoid)\",6006,9122,7046,8901,0,0,,,,,,,,"));
			coordCSV.Add(4007, GetCoordCSV("4007,Unknown datum based upon the Clarke 1858 ellipsoid,6007,\"Not specified (based on Clarke 1858 ellipsoid)\",6007,9122,7007,8901,0,0,,,,,,,,"));
			coordCSV.Add(4008, GetCoordCSV("4008,Unknown datum based upon the Clarke 1866 ellipsoid,6008,\"Not specified (based on Clarke 1866 ellipsoid)\",6008,9122,7008,8901,0,0,,,,,,,,"));
			coordCSV.Add(4009, GetCoordCSV("4009,Unknown datum based upon the Clarke 1866 Michigan ellipsoid,6009,\"Not specified (based on Clarke 1866 Michigan ellipsoid)\",6009,9122,7009,8901,0,0,,,,,,,,"));
			coordCSV.Add(4010, GetCoordCSV("4010,\"Unknown datum based upon the Clarke 1880 (Benoit) ellipsoid\",6010,\"Not specified (based on Clarke 1880 (Benoit) ellipsoid)\",6010,9122,7010,8901,0,0,,,,,,,,"));
			coordCSV.Add(4011, GetCoordCSV("4011,\"Unknown datum based upon the Clarke 1880 (IGN) ellipsoid\",6011,\"Not specified (based on Clarke 1880 (IGN) ellipsoid)\",6011,9122,7011,8901,0,0,,,,,,,,"));
			coordCSV.Add(4012, GetCoordCSV("4012,\"Unknown datum based upon the Clarke 1880 (RGS) ellipsoid\",6012,\"Not specified (based on Clarke 1880 (RGS) ellipsoid)\",6012,9122,7012,8901,0,0,,,,,,,,"));
			coordCSV.Add(4013, GetCoordCSV("4013,\"Unknown datum based upon the Clarke 1880 (Arc) ellipsoid\",6013,\"Not specified (based on Clarke 1880 (Arc) ellipsoid)\",6013,9122,7013,8901,0,0,,,,,,,,"));
			coordCSV.Add(4014, GetCoordCSV("4014,\"Unknown datum based upon the Clarke 1880 (SGA 1922) ellipsoid\",6014,\"Not specified (based on Clarke 1880 (SGA 1922) ellipsoid)\",6014,9122,7014,8901,0,0,,,,,,,,"));
			coordCSV.Add(4015, GetCoordCSV("4015,\"Unknown datum based upon the Everest 1830 (1937 Adjustment) ellipsoid\",6015,\"Not specified (based on Everest 1830 (1937 Adjustment) ellipsoid)\",6015,9122,7015,8901,0,0,,,,,,,,"));
			coordCSV.Add(4016, GetCoordCSV("4016,\"Unknown datum based upon the Everest 1830 (1967 Definition) ellipsoid\",6016,\"Not specified (based on Everest 1830 (1967 Definition) ellipsoid)\",6016,9122,7016,8901,0,0,,,,,,,,"));
			coordCSV.Add(4018, GetCoordCSV("4018,Unknown datum based upon the Everest 1830 Modified ellipsoid,6018,\"Not specified (based on Everest 1830 Modified ellipsoid)\",6018,9122,7018,8901,0,0,,,,,,,,"));
			coordCSV.Add(4019, GetCoordCSV("4019,Unknown datum based upon the GRS 1980 ellipsoid,6019,\"Not specified (based on GRS 1980 ellipsoid)\",6019,9122,7019,8901,0,0,,,,,,,,"));
			coordCSV.Add(4020, GetCoordCSV("4020,Unknown datum based upon the Helmert 1906 ellipsoid,6020,\"Not specified (based on Helmert 1906 ellipsoid)\",6020,9122,7020,8901,0,0,,,,,,,,"));
			coordCSV.Add(4021, GetCoordCSV("4021,Unknown datum based upon the Indonesian National Spheroid,6021,\"Not specified (based on Indonesian National Spheroid)\",6021,9122,7021,8901,0,0,,,,,,,,"));
			coordCSV.Add(4022, GetCoordCSV("4022,Unknown datum based upon the International 1924 ellipsoid,6022,\"Not specified (based on International 1924 ellipsoid)\",6022,9122,7022,8901,0,0,,,,,,,,"));
			coordCSV.Add(4024, GetCoordCSV("4024,Unknown datum based upon the Krassowsky 1940 ellipsoid,6024,\"Not specified (based on Krassowsky 1940 ellipsoid)\",6024,9122,7024,8901,0,0,,,,,,,,"));
			coordCSV.Add(4025, GetCoordCSV("4025,Unknown datum based upon the NWL 9D ellipsoid,6025,\"Not specified (based on NWL 9D ellipsoid)\",6025,9122,7025,8901,0,0,,,,,,,,"));
			coordCSV.Add(4027, GetCoordCSV("4027,Unknown datum based upon the Plessis 1817 ellipsoid,6027,\"Not specified (based on Plessis 1817 ellipsoid)\",6027,9122,7027,8901,0,0,,,,,,,,"));
			coordCSV.Add(4028, GetCoordCSV("4028,Unknown datum based upon the Struve 1860 ellipsoid,6028,\"Not specified (based on Struve 1860 ellipsoid)\",6028,9122,7028,8901,0,0,,,,,,,,"));
			coordCSV.Add(4029, GetCoordCSV("4029,Unknown datum based upon the War Office ellipsoid,6029,\"Not specified (based on War Office ellipsoid)\",6029,9122,7029,8901,0,0,,,,,,,,"));
			coordCSV.Add(4030, GetCoordCSV("4030,Unknown datum based upon the WGS 84 ellipsoid,6030,\"Not specified (based on WGS 84 ellipsoid)\",6030,9122,7030,8901,0,0,,,,,,,,"));
			coordCSV.Add(4031, GetCoordCSV("4031,Unknown datum based upon the GEM 10C ellipsoid,6031,\"Not specified (based on GEM 10C ellipsoid)\",6031,9122,7031,8901,0,0,,,,,,,,"));
			coordCSV.Add(4032, GetCoordCSV("4032,Unknown datum based upon the OSU86F ellipsoid,6032,\"Not specified (based on OSU86F ellipsoid)\",6032,9122,7032,8901,0,0,,,,,,,,"));
			coordCSV.Add(4033, GetCoordCSV("4033,Unknown datum based upon the OSU91A ellipsoid,6033,\"Not specified (based on OSU91A ellipsoid)\",6033,9122,7033,8901,0,0,,,,,,,,"));
			coordCSV.Add(4034, GetCoordCSV("4034,Unknown datum based upon the Clarke 1880 ellipsoid,6034,\"Not specified (based on Clarke 1880 ellipsoid)\",6034,9122,7034,8901,0,0,,,,,,,,"));
			coordCSV.Add(4035, GetCoordCSV("4035,Unknown datum based upon the Authalic Sphere,6035,\"Not specified (based on Authalic Sphere)\",6035,9108,7035,8901,0,1,,,,,,,,"));
			coordCSV.Add(4036, GetCoordCSV("4036,Unknown datum based upon the GRS 1967 ellipsoid,6036,\"Not specified (based on GRS 1967 ellipsoid)\",6036,9122,7036,8901,0,0,,,,,,,,"));
			coordCSV.Add(4041, GetCoordCSV("4041,Unknown datum based upon the Average Terrestrial System 1977 ellipsoid,6041,\"Not specified (based on Average Terrestrial System 1977 ellipsoid)\",6041,9122,7041,8901,0,0,,,,,,,,"));
			coordCSV.Add(4042, GetCoordCSV("4042,\"Unknown datum based upon the Everest (1830 Definition) ellipsoid\",6042,\"Not specified (based on Everest (1830 Definition) ellipsoid)\",6042,9122,7042,8901,0,0,,,,,,,,"));
			coordCSV.Add(4043, GetCoordCSV("4043,Unknown datum based upon the WGS 72 ellipsoid,6043,\"Not specified (based on WGS 72 ellipsoid)\",6043,9122,7043,8901,0,0,,,,,,,,"));
			coordCSV.Add(4044, GetCoordCSV("4044,\"Unknown datum based upon the Everest 1830 (1962 Definition) ellipsoid\",6044,\"Not specified (based on Everest 1830 (1962 Definition) ellipsoid)\",6044,9122,7044,8901,0,0,,,,,,,,"));
			coordCSV.Add(4045, GetCoordCSV("4045,\"Unknown datum based upon the Everest 1830 (1975 Definition) ellipsoid\",6045,\"Not specified (based on Everest 1830 (1975 Definition) ellipsoid)\",6045,9122,7045,8901,0,0,,,,,,,,"));
			coordCSV.Add(4047, GetCoordCSV("4047,Unspecified datum based upon the GRS 1980 Authalic Sphere,6047,\"Not specified (based on GRS 1980 Authalic Sphere)\",6047,9122,7048,8901,0,0,,,,,,,,"));
			coordCSV.Add(4052, GetCoordCSV("4052,Unspecified datum based upon the Clarke 1866 Authalic Sphere,6052,\"Not specified (based on Clarke 1866 Authalic Sphere)\",6052,9122,7052,8901,0,0,,,,,,,,"));
			coordCSV.Add(4053, GetCoordCSV("4053,Unspecified datum based upon the International 1924 Authalic Sphere,6053,\"Not specified (based on International 1924 Authalic Sphere)\",6053,9122,7057,8901,0,0,,,,,,,,"));
			coordCSV.Add(4054, GetCoordCSV("4054,Unspecified datum based upon the Hughes 1980 ellipsoid,6054,\"Not specified (based on Hughes 1980 ellipsoid)\",6054,9122,7058,8901,0,0,,,,,,,,"));
			coordCSV.Add(4120, GetCoordCSV("4120,Greek,6120,Greek,6120,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4121, GetCoordCSV("4121,GGRS87,6121,Greek Geodetic Reference System 1987,6121,9122,7019,8901,1,0,9603,-199.87,74.79,246.62,,,,"));
			coordCSV.Add(4122, GetCoordCSV("4122,ATS77,6122,Average Terrestrial System 1977,6122,9122,7041,8901,1,0,,,,,,,,"));
			//coordCSV.Add(4123, GetCoordCSV("4123,KKJ,6123,\"Kartastokoordinaattijarjestelma (1966)\",6123,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4124, GetCoordCSV("4124,RT90,6124,Rikets koordinatsystem 1990,6124,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4125, GetCoordCSV("4125,Samboja,6125,Samboja,6125,9108,7004,8901,1,1,9603,-404.78,685.68,45.47,,,,"));
			coordCSV.Add(4126, GetCoordCSV("4126,\"LKS94 (ETRS89)\",6126,\"Lithuania 1994 (ETRS89)\",6126,9108,7019,8901,1,1,,,,,,,,"));
			coordCSV.Add(4127, GetCoordCSV("4127,Tete,6127,Tete,6127,9122,7008,8901,1,0,,,,,,,,"));
			coordCSV.Add(4128, GetCoordCSV("4128,Madzansua,6128,Madzansua,6128,9122,7008,8901,1,0,,,,,,,,"));
			coordCSV.Add(4129, GetCoordCSV("4129,Observatario,6129,Observatario,6129,9122,7008,8901,1,0,,,,,,,,"));
			coordCSV.Add(4130, GetCoordCSV("4130,Moznet,6130,\"Moznet (ITRF94)\",6130,9122,7030,8901,1,0,9607,0,0,0,0,0,0,0"));
			coordCSV.Add(4131, GetCoordCSV("4131,Indian 1960,6131,Indian 1960,6131,9122,7015,8901,1,0,,,,,,,,"));
			coordCSV.Add(4132, GetCoordCSV("4132,FD58,6132,Final Datum 1958,6132,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4133, GetCoordCSV("4133,EST92,6133,Estonia 1992,6133,9122,7019,8901,1,0,9607,0.055,-0.541,-0.185,-0.0183,0.0003,0.007,-0.014"));
			coordCSV.Add(4134, GetCoordCSV("4134,PDO Survey Datum 1993,6134,PDO Survey Datum 1993,6134,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4135, GetCoordCSV("4135,Old Hawaiian,6135,Old Hawaiian,6135,9122,7008,8901,1,0,,,,,,,,"));
			coordCSV.Add(4136, GetCoordCSV("4136,St. Lawrence Island,6136,St. Lawrence Island,6136,9122,7008,8901,1,0,,,,,,,,"));
			coordCSV.Add(4137, GetCoordCSV("4137,St. Paul Island,6137,St. Paul Island,6137,9122,7008,8901,1,0,,,,,,,,"));
			coordCSV.Add(4138, GetCoordCSV("4138,St. George Island,6138,St. George Island,6138,9122,7008,8901,1,0,,,,,,,,"));
			coordCSV.Add(4139, GetCoordCSV("4139,Puerto Rico,6139,Puerto Rico,6139,9122,7008,8901,1,0,9603,11,72,-101,,,,"));
			coordCSV.Add(4140, GetCoordCSV("4140,\"NAD83(CSRS98)\",6140,NAD83 Canadian Spatial Reference System,6140,9108,7019,8901,1,1,9603,0,0,0,,,,"));
			coordCSV.Add(4141, GetCoordCSV("4141,Israel,6141,Israel,6141,9122,7019,8901,1,0,9603,-48,55,52,,,,"));
			coordCSV.Add(4142, GetCoordCSV("4142,Locodjo 1965,6142,Locodjo 1965,6142,9122,7012,8901,1,0,9603,-125,53,467,,,,"));
			coordCSV.Add(4143, GetCoordCSV("4143,Abidjan 1987,6143,Abidjan 1987,6143,9122,7012,8901,1,0,9603,-124.76,53,466.79,,,,"));
			coordCSV.Add(4144, GetCoordCSV("4144,Kalianpur 1937,6144,Kalianpur 1937,6144,9122,7015,8901,1,0,,,,,,,,"));
			coordCSV.Add(4145, GetCoordCSV("4145,Kalianpur 1962,6145,Kalianpur 1962,6145,9122,7044,8901,1,0,,,,,,,,"));
			coordCSV.Add(4146, GetCoordCSV("4146,Kalianpur 1975,6146,Kalianpur 1975,6146,9122,7045,8901,1,0,9603,295,736,257,,,,"));
			coordCSV.Add(4147, GetCoordCSV("4147,Hanoi 1972,6147,Hanoi 1972,6147,9122,7024,8901,1,0,9603,-17.51,-108.32,-62.39,,,,"));
			coordCSV.Add(4148, GetCoordCSV("4148,Hartebeesthoek94,6148,Hartebeesthoek94,6148,9122,7030,8901,1,0,9603,0,0,0,,,,"));
			//coordCSV.Add(4149, GetCoordCSV("4149,CH1903,6149,CH1903,6149,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4150, GetCoordCSV("4150,\"CH1903+\",6150,\"CH1903+\",6150,9122,7004,8901,1,0,9603,674.374,15.056,405.346,,,,"));
			coordCSV.Add(4151, GetCoordCSV("4151,CHTRF95,6151,Swiss Terrestrial Reference Frame 1995,6151,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4152, GetCoordCSV("4152,\"NAD83(HARN)\",6152,\"NAD83 (High Accuracy Regional Network)\",6152,9122,7019,8901,1,0,,,,,,,,"));
			coordCSV.Add(4153, GetCoordCSV("4153,Rassadiran,6153,Rassadiran,6153,9122,7022,8901,1,0,9603,-133.63,-157.5,-158.62,,,,"));
			coordCSV.Add(4154, GetCoordCSV("4154,\"ED50(ED77)\",6154,\"European Datum 1950(1977)\",6154,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4155, GetCoordCSV("4155,Dabola 1981,6155,Dabola 1981,6155,9122,7011,8901,1,0,9603,-83,37,124,,,,"));
			coordCSV.Add(4156, GetCoordCSV("4156,S-JTSK,6156,Jednotne Trigonometricke Site Katastralni,6156,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4157, GetCoordCSV("4157,Mount Dillon,6157,Mount Dillon,6157,9122,7007,8901,1,0,,,,,,,,"));
			coordCSV.Add(4158, GetCoordCSV("4158,Naparima 1955,6158,Naparima 1955,6158,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4159, GetCoordCSV("4159,ELD79,6159,European Libyan Datum 1979,6159,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4160, GetCoordCSV("4160,Chos Malal 1914,6160,Chos Malal 1914,6160,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4161, GetCoordCSV("4161,Pampa del Castillo,6161,Pampa del Castillo,6161,9122,7022,8901,1,0,9603,27.5,14,186.4,,,,"));
			coordCSV.Add(4162, GetCoordCSV("4162,Korean 1985,6162,Korean Datum 1985,6162,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4163, GetCoordCSV("4163,Yemen NGN96,6163,Yemen National Geodetic Network 1996,6163,9122,7030,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4164, GetCoordCSV("4164,South Yemen,6164,South Yemen,6164,9122,7024,8901,1,0,9603,-76,-138,67,,,,"));
			coordCSV.Add(4165, GetCoordCSV("4165,Bissau,6165,Bissau,6165,9122,7022,8901,1,0,9603,-173,253,27,,,,"));
			coordCSV.Add(4166, GetCoordCSV("4166,Korean 1995,6166,Korean Datum 1995,6166,9122,7030,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4167, GetCoordCSV("4167,NZGD2000,6167,New Zealand Geodetic Datum 2000,6167,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4168, GetCoordCSV("4168,Accra,6168,Accra,6168,9122,7029,8901,1,0,9603,-199,32,322,,,,"));
			coordCSV.Add(4169, GetCoordCSV("4169,American Samoa 1962,6169,American Samoa 1962,6169,9122,7008,8901,1,0,9603,-115,118,426,,,,"));
			coordCSV.Add(4170, GetCoordCSV("4170,SIRGAS,6170,Sistema de Referencia Geocentrico para America del Sur 1995,6170,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4171, GetCoordCSV("4171,RGF93,6171,Reseau Geodesique Francais 1993,6171,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4172, GetCoordCSV("4172,POSGAR,6172,Posiciones Geodesicas Argentinas,6172,9108,7019,8901,1,1,9603,0,0,0,,,,"));
			coordCSV.Add(4173, GetCoordCSV("4173,IRENET95,6173,IRENET95,6173,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4174, GetCoordCSV("4174,Sierra Leone 1924,6174,Sierra Leone Colony 1924,6174,9122,7029,8901,1,0,,,,,,,,"));
			coordCSV.Add(4175, GetCoordCSV("4175,Sierra Leone 1968,6175,Sierra Leone 1968,6175,9122,7012,8901,1,0,9603,-88,4,101,,,,"));
			coordCSV.Add(4176, GetCoordCSV("4176,Australian Antarctic,6176,Australian Antarctic Datum 1998,6176,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4178, GetCoordCSV("4178,\"Pulkovo 1942(83)\",6178,\"Pulkovo 1942/83\",6178,9122,7024,8901,1,0,9607,24,-123,-94,-0.02,0.25,0.13,1.1"));
			coordCSV.Add(4179, GetCoordCSV("4179,\"Pulkovo 1942(58)\",6179,\"Pulkovo 1942/58\",6179,9122,7024,8901,1,0,9606,33.4,-146.6,-76.3,-0.359,-0.053,0.844,-0.84"));
			coordCSV.Add(4180, GetCoordCSV("4180,EST97,6180,Estonia 1997,6180,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4181, GetCoordCSV("4181,Luxembourg 1930,6181,Luxembourg 1930,6181,9122,7022,8901,1,0,9606,-193,13.7,-39.3,-0.41,-2.933,2.688,0.43"));
			coordCSV.Add(4182, GetCoordCSV("4182,Azores Occidental 1939,6182,Azores Occidental Islands 1939,6182,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4183, GetCoordCSV("4183,Azores Central 1948,6183,Azores Central Islands 1948,6183,9122,7022,8901,1,0,9603,-104,167,-38,,,,"));
			coordCSV.Add(4184, GetCoordCSV("4184,Azores Oriental 1940,6184,Azores Oriental Islands 1940,6184,9122,7022,8901,1,0,9603,-203,141,53,,,,"));
			coordCSV.Add(4185, GetCoordCSV("4185,Madeira 1936,6185,Madeira 1936,6185,9108,7022,8901,1,1,,,,,,,,"));
			coordCSV.Add(4188, GetCoordCSV("4188,OSNI 1952,6188,OSNI 1952,6188,9122,7001,8901,1,0,9606,482.5,-130.6,564.6,-1.042,-0.214,-0.631,8.15"));
			coordCSV.Add(4189, GetCoordCSV("4189,REGVEN,6189,Red Geodesica Venezolana,6189,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4190, GetCoordCSV("4190,POSGAR 98,6190,Posiciones Geodesicas Argentinas 1998,6190,9122,7019,8901,1,0,,,,,,,,"));
			coordCSV.Add(4191, GetCoordCSV("4191,Albanian 1987,6191,Albanian 1987,6191,9122,7024,8901,1,0,,,,,,,,"));
			coordCSV.Add(4192, GetCoordCSV("4192,Douala 1948,6192,Douala 1948,6192,9122,7022,8901,1,0,9603,-206.1,-174.7,-87.7,,,,"));
			coordCSV.Add(4193, GetCoordCSV("4193,Manoca 1962,6193,Manoca 1962,6193,9122,7011,8901,1,0,9603,-70.9,-151.8,-41.4,,,,"));
			coordCSV.Add(4194, GetCoordCSV("4194,Qornoq 1927,6194,Qornoq 1927,6194,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4195, GetCoordCSV("4195,Scoresbysund 1952,6195,Scoresbysund 1952,6195,9122,7022,8901,1,0,9606,105,326,-102.5,0,0,0.814,-0.6"));
			coordCSV.Add(4196, GetCoordCSV("4196,Ammassalik 1958,6196,Ammassalik 1958,6196,9122,7022,8901,1,0,9606,-45,417,-3.5,0,0,0.814,-0.6"));
			coordCSV.Add(4197, GetCoordCSV("4197,Garoua,6197,Garoua,6197,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4198, GetCoordCSV("4198,Kousseri,6198,Kousseri,6198,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4199, GetCoordCSV("4199,Egypt 1930,6199,Egypt 1930,6199,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4200, GetCoordCSV("4200,Pulkovo 1995,6200,Pulkovo 1995,6200,9122,7024,8901,1,0,,,,,,,,"));
			coordCSV.Add(4201, GetCoordCSV("4201,Adindan,6201,Adindan,6201,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4202, GetCoordCSV("4202,AGD66,6202,Australian Geodetic Datum 1966,6202,9122,7003,8901,1,0,,,,,,,,"));
			coordCSV.Add(4203, GetCoordCSV("4203,AGD84,6203,Australian Geodetic Datum 1984,6203,9122,7003,8901,1,0,,,,,,,,"));
			coordCSV.Add(4204, GetCoordCSV("4204,Ain el Abd,6204,Ain el Abd 1970,6204,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4205, GetCoordCSV("4205,Afgooye,6205,Afgooye,6205,9122,7024,8901,1,0,9603,-43,-163,45,,,,"));
			coordCSV.Add(4206, GetCoordCSV("4206,Agadez,6206,Agadez,6206,9122,7011,8901,1,0,,,,,,,,"));
			coordCSV.Add(4207, GetCoordCSV("4207,Lisbon,6207,Lisbon 1937,6207,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4208, GetCoordCSV("4208,Aratu,6208,Aratu,6208,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4209, GetCoordCSV("4209,Arc 1950,6209,Arc 1950,6209,9122,7013,8901,1,0,,,,,,,,"));
			coordCSV.Add(4210, GetCoordCSV("4210,Arc 1960,6210,Arc 1960,6210,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4211, GetCoordCSV("4211,Batavia,6211,Batavia,6211,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4212, GetCoordCSV("4212,Barbados 1938,6212,Barbados 1938,6212,9122,7012,8901,1,0,9603,31.95,300.99,419.19,,,,"));
			coordCSV.Add(4213, GetCoordCSV("4213,Beduaram,6213,Beduaram,6213,9122,7011,8901,1,0,9603,-106,-87,188,,,,"));
			coordCSV.Add(4214, GetCoordCSV("4214,Beijing 1954,6214,Beijing 1954,6214,9122,7024,8901,1,0,,,,,,,,"));
			coordCSV.Add(4215, GetCoordCSV("4215,Belge 1950,6215,Reseau National Belge 1950,6215,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4216, GetCoordCSV("4216,Bermuda 1957,6216,Bermuda 1957,6216,9122,7008,8901,1,0,9603,-73,213,296,,,,"));
			//coordCSV.Add(4218, GetCoordCSV("4218,Bogota 1975,6218,Bogota 1975,6218,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4219, GetCoordCSV("4219,Bukit Rimpah,6219,Bukit Rimpah,6219,9122,7004,8901,1,0,9603,-384,664,-48,,,,"));
			coordCSV.Add(4220, GetCoordCSV("4220,Camacupa,6220,Camacupa,6220,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4221, GetCoordCSV("4221,Campo Inchauspe,6221,Campo Inchauspe,6221,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4222, GetCoordCSV("4222,Cape,6222,Cape,6222,9122,7013,8901,1,0,,,,,,,,"));
			coordCSV.Add(4223, GetCoordCSV("4223,Carthage,6223,Carthage,6223,9122,7011,8901,1,0,,,,,,,,"));
			coordCSV.Add(4224, GetCoordCSV("4224,Chua,6224,Chua,6224,9122,7022,8901,1,0,9603,-134,229,-29,,,,"));
			coordCSV.Add(4225, GetCoordCSV("4225,Corrego Alegre,6225,Corrego Alegre,6225,9122,7022,8901,1,0,9603,-206,172,-6,,,,"));
			coordCSV.Add(4226, GetCoordCSV("4226,\"Cote d'Ivoire\",6226,\"Cote d'Ivoire\",6226,9108,7011,8901,1,1,,,,,,,,"));
			coordCSV.Add(4227, GetCoordCSV("4227,Deir ez Zor,6227,Deir ez Zor,6227,9122,7011,8901,1,0,,,,,,,,"));
			coordCSV.Add(4228, GetCoordCSV("4228,Douala,6228,Douala,6228,9108,7011,8901,1,1,,,,,,,,"));
			coordCSV.Add(4229, GetCoordCSV("4229,Egypt 1907,6229,Egypt 1907,6229,9122,7020,8901,1,0,,,,,,,,"));
			coordCSV.Add(4230, GetCoordCSV("4230,ED50,6230,European Datum 1950,6230,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4231, GetCoordCSV("4231,ED87,6231,European Datum 1987,6231,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4232, GetCoordCSV("4232,Fahud,6232,Fahud,6232,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4233, GetCoordCSV("4233,Gandajika 1970,6233,Gandajika 1970,6233,9122,7022,8901,1,1,9603,-133,-321,50,,,,"));
			coordCSV.Add(4234, GetCoordCSV("4234,Garoua,6234,Garoua,6234,9108,7011,8901,1,1,,,,,,,,"));
			coordCSV.Add(4235, GetCoordCSV("4235,Guyane Francaise,6235,Guyane Francaise,6235,9108,7022,8901,1,1,,,,,,,,"));
			coordCSV.Add(4236, GetCoordCSV("4236,Hu Tzu Shan,6236,Hu Tzu Shan,6236,9122,7022,8901,1,0,9603,-637,-549,-203,,,,"));
			coordCSV.Add(4237, GetCoordCSV("4237,HD72,6237,Hungarian Datum 1972,6237,9122,7036,8901,1,0,,,,,,,,"));
			coordCSV.Add(4238, GetCoordCSV("4238,ID74,6238,Indonesian Datum 1974,6238,9122,7021,8901,1,0,,,,,,,,"));
			coordCSV.Add(4239, GetCoordCSV("4239,Indian 1954,6239,Indian 1954,6239,9122,7015,8901,1,0,9603,217,823,299,,,,"));
			coordCSV.Add(4240, GetCoordCSV("4240,Indian 1975,6240,Indian 1975,6240,9122,7015,8901,1,0,,,,,,,,"));
			coordCSV.Add(4241, GetCoordCSV("4241,Jamaica 1875,6241,Jamaica 1875,6241,9122,7034,8901,1,0,,,,,,,,"));
			coordCSV.Add(4242, GetCoordCSV("4242,JAD69,6242,Jamaica 1969,6242,9122,7008,8901,1,0,,,,,,,,"));
			coordCSV.Add(4243, GetCoordCSV("4243,Kalianpur 1880,6243,Kalianpur 1880,6243,9122,7042,8901,1,0,,,,,,,,"));
			coordCSV.Add(4244, GetCoordCSV("4244,Kandawala,6244,Kandawala,6244,9122,7015,8901,1,0,9603,-97,787,86,,,,"));
			coordCSV.Add(4245, GetCoordCSV("4245,Kertau 1968,6245,Kertau 1968,6245,9122,7018,8901,1,0,9603,-11,851,5,,,,"));
			coordCSV.Add(4246, GetCoordCSV("4246,KOC,6246,Kuwait Oil Company,6246,9122,7012,8901,1,0,9603,-294.7,-200.1,525.5,,,,"));
			coordCSV.Add(4247, GetCoordCSV("4247,La Canoa,6247,La Canoa,6247,9122,7022,8901,1,0,9603,-273.5,110.6,-357.9,,,,"));
			coordCSV.Add(4248, GetCoordCSV("4248,PSAD56,6248,Provisional South American Datum 1956,6248,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4249, GetCoordCSV("4249,Lake,6249,Lake,6249,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4250, GetCoordCSV("4250,Leigon,6250,Leigon,6250,9122,7012,8901,1,0,9603,-130,29,364,,,,"));
			coordCSV.Add(4251, GetCoordCSV("4251,Liberia 1964,6251,Liberia 1964,6251,9122,7012,8901,1,0,9603,-90,40,88,,,,"));
			coordCSV.Add(4252, GetCoordCSV("4252,Lome,6252,Lome,6252,9122,7011,8901,1,0,,,,,,,,"));
			coordCSV.Add(4253, GetCoordCSV("4253,Luzon 1911,6253,Luzon 1911,6253,9122,7008,8901,1,0,,,,,,,,"));
			coordCSV.Add(4254, GetCoordCSV("4254,Hito XVIII 1963,6254,Hito XVIII 1963,6254,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4255, GetCoordCSV("4255,Herat North,6255,Herat North,6255,9122,7022,8901,1,0,9603,-333,-222,114,,,,"));
			coordCSV.Add(4256, GetCoordCSV("4256,Mahe 1971,6256,Mahe 1971,6256,9122,7012,8901,1,0,9603,41,-220,-134,,,,"));
			coordCSV.Add(4257, GetCoordCSV("4257,Makassar,6257,Makassar,6257,9122,7004,8901,1,0,9603,-587.8,519.75,145.76,,,,"));
			coordCSV.Add(4258, GetCoordCSV("4258,ETRS89,6258,European Terrestrial Reference System 1989,6258,9122,7019,8901,1,0,,,,,,,,"));
			coordCSV.Add(4259, GetCoordCSV("4259,Malongo 1987,6259,Malongo 1987,6259,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4260, GetCoordCSV("4260,Manoca,6260,Manoca,6260,9108,7012,8901,1,1,9603,-70.9,-151.8,-41.4,,,,"));
			coordCSV.Add(4261, GetCoordCSV("4261,Merchich,6261,Merchich,6261,9122,7011,8901,1,0,9603,31,146,47,,,,"));
			coordCSV.Add(4262, GetCoordCSV("4262,Massawa,6262,Massawa,6262,9122,7004,8901,1,0,9603,639,405,60,,,,"));
			coordCSV.Add(4263, GetCoordCSV("4263,Minna,6263,Minna,6263,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4264, GetCoordCSV("4264,Mhast,6264,Mhast,6264,9122,7022,8901,1,1,9603,-252.95,-4.11,-96.38,,,,"));
			coordCSV.Add(4265, GetCoordCSV("4265,Monte Mario,6265,Monte Mario,6265,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4266, GetCoordCSV("4266,\"M'poraloko\",6266,\"M'poraloko\",6266,9122,7011,8901,1,0,,,,,,,,"));
			coordCSV.Add(4267, GetCoordCSV("4267,NAD27,6267,North American Datum 1927,6267,9122,7008,8901,1,0,,,,,,,,"));
			coordCSV.Add(4268, GetCoordCSV("4268,NAD27 Michigan,6268,NAD Michigan,6268,9122,7009,8901,1,0,,,,,,,,"));
			coordCSV.Add(4269, GetCoordCSV("4269,NAD83,6269,North American Datum 1983,6269,9122,7019,8901,1,0,,,,,,,,"));
			coordCSV.Add(4270, GetCoordCSV("4270,Nahrwan 1967,6270,Nahrwan 1967,6270,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4271, GetCoordCSV("4271,Naparima 1972,6271,Naparima 1972,6271,9122,7022,8901,1,0,,,,,,,,"));
			//coordCSV.Add(4272, GetCoordCSV("4272,NZGD49,6272,New Zealand Geodetic Datum 1949,6272,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4273, GetCoordCSV("4273,NGO 1948,6273,NGO 1948,6273,9122,7005,8901,1,0,9606,278.3,93,474.5,7.889,0.05,-6.61,6.21"));
			coordCSV.Add(4274, GetCoordCSV("4274,Datum 73,6274,Datum 73,6274,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4275, GetCoordCSV("4275,NTF,6275,Nouvelle Triangulation Francaise,6275,9122,7011,8901,1,0,9603,-168,-60,320,,,,"));
			coordCSV.Add(4276, GetCoordCSV("4276,NSWC 9Z-2,6276,NSWC 9Z-2,6276,9122,7025,8901,1,0,,,,,,,,"));
			coordCSV.Add(4277, GetCoordCSV("4277,OSGB 1936,6277,OSGB 1936,6277,9122,7001,8901,1,0,,,,,,,,"));
			coordCSV.Add(4278, GetCoordCSV("4278,OSGB70,6278,\"OSGB 1970 (SN)\",6278,9122,7001,8901,1,0,,,,,,,,"));
			coordCSV.Add(4279, GetCoordCSV("4279,\"OS(SN)80\",6279,\"OS (SN) 1980\",6279,9122,7001,8901,1,0,,,,,,,,"));
			coordCSV.Add(4280, GetCoordCSV("4280,Padang,6280,Padang 1884,6280,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4281, GetCoordCSV("4281,Palestine 1923,6281,Palestine 1923,6281,9122,7010,8901,1,0,9606,-275.7224,94.7824,340.8944,-8.001,-4.42,-11.821,1"));
			coordCSV.Add(4282, GetCoordCSV("4282,Pointe Noire,6282,Congo 1960 Pointe Noire,6282,9122,7011,8901,1,0,,,,,,,,"));
			coordCSV.Add(4283, GetCoordCSV("4283,GDA94,6283,Geocentric Datum of Australia 1994,6283,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4284, GetCoordCSV("4284,Pulkovo 1942,6284,Pulkovo 1942,6284,9122,7024,8901,1,0,,,,,,,,"));
			coordCSV.Add(4285, GetCoordCSV("4285,Qatar 1974,6285,Qatar 1974,6285,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4286, GetCoordCSV("4286,Qatar 1948,6286,Qatar 1948,6286,9122,7020,8901,1,0,,,,,,,,"));
			coordCSV.Add(4287, GetCoordCSV("4287,Qornoq,6287,Qornoq,6287,9108,7022,8901,1,1,9603,164,138,-189,,,,"));
			coordCSV.Add(4288, GetCoordCSV("4288,Loma Quintana,6288,Loma Quintana,6288,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4289, GetCoordCSV("4289,Amersfoort,6289,Amersfoort,6289,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4291, GetCoordCSV("4291,SAD69,6291,South American Datum 1969,6291,9108,7036,8901,1,1,,,,,,,,"));
			coordCSV.Add(4292, GetCoordCSV("4292,Sapper Hill 1943,6292,Sapper Hill 1943,6292,9122,7022,8901,1,0,9603,-355,21,72,,,,"));
			coordCSV.Add(4293, GetCoordCSV("4293,Schwarzeck,6293,Schwarzeck,6293,9122,7046,8901,1,0,,,,,,,,"));
			coordCSV.Add(4294, GetCoordCSV("4294,Segora,6294,Segora,6294,9108,7004,8901,1,1,,,,,,,,"));
			coordCSV.Add(4295, GetCoordCSV("4295,Serindung,6295,Serindung,6295,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4296, GetCoordCSV("4296,Sudan,6296,Sudan,6296,9108,7011,8901,1,1,,,,,,,,"));
			coordCSV.Add(4297, GetCoordCSV("4297,Tananarive,6297,Tananarive 1925,6297,9122,7022,8901,1,0,9603,-189,-242,-91,,,,"));
			coordCSV.Add(4298, GetCoordCSV("4298,Timbalai 1948,6298,Timbalai 1948,6298,9122,7016,8901,1,0,,,,,,,,"));
			coordCSV.Add(4299, GetCoordCSV("4299,TM65,6299,TM65,6299,9122,7002,8901,1,0,,,,,,,,"));
			coordCSV.Add(4300, GetCoordCSV("4300,TM75,6300,Geodetic Datum of 1965,6300,9122,7002,8901,1,0,,,,,,,,"));
			coordCSV.Add(4301, GetCoordCSV("4301,Tokyo,6301,Tokyo,6301,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4302, GetCoordCSV("4302,Trinidad 1903,6302,Trinidad 1903,6302,9122,7007,8901,1,0,,,,,,,,"));
			coordCSV.Add(4303, GetCoordCSV("4303,\"TC(1948)\",6303,Trucial Coast 1948,6303,9122,7020,8901,1,0,,,,,,,,"));
			coordCSV.Add(4304, GetCoordCSV("4304,Voirol 1875,6304,Voirol 1875,6304,9122,7011,8901,1,0,9603,-73,-247,227,,,,"));
			coordCSV.Add(4306, GetCoordCSV("4306,Bern 1938,6306,Bern 1938,6306,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4307, GetCoordCSV("4307,Nord Sahara 1959,6307,Nord Sahara 1959,6307,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4308, GetCoordCSV("4308,RT38,6308,Stockholm 1938,6308,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4309, GetCoordCSV("4309,Yacare,6309,Yacare,6309,9122,7022,8901,1,0,9603,-155,171,37,,,,"));
			coordCSV.Add(4310, GetCoordCSV("4310,Yoff,6310,Yoff,6310,9122,7011,8901,1,0,,,,,,,,"));
			coordCSV.Add(4311, GetCoordCSV("4311,Zanderij,6311,Zanderij,6311,9122,7022,8901,1,0,9603,-265,120,-358,,,,"));
			//coordCSV.Add(4312, GetCoordCSV("4312,MGI,6312,Militar-Geographische Institut,6312,9122,7004,8901,1,0,,,,,,,,"));
			//coordCSV.Add(4313, GetCoordCSV("4313,Belge 1972,6313,Reseau National Belge 1972,6313,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4314, GetCoordCSV("4314,DHDN,6314,Deutsches Hauptdreiecksnetz,6314,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4315, GetCoordCSV("4315,Conakry 1905,6315,Conakry 1905,6315,9122,7011,8901,1,0,9603,-23,259,-9,,,,"));
			coordCSV.Add(4316, GetCoordCSV("4316,Dealul Piscului 1933,6316,Dealul Piscului 1933,6316,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4317, GetCoordCSV("4317,Dealul Piscului 1970,6317,Dealul Piscului 1970,6317,9122,7024,8901,1,0,,,,,,,,"));
			coordCSV.Add(4318, GetCoordCSV("4318,NGN,6318,National Geodetic Network,6318,9122,7030,8901,1,0,9603,-3.2,-5.7,2.8,,,,"));
			coordCSV.Add(4319, GetCoordCSV("4319,KUDAMS,6319,Kuwait Utility,6319,9122,7019,8901,1,0,,,,,,,,"));
			coordCSV.Add(4322, GetCoordCSV("4322,WGS 72,6322,World Geodetic System 1972,6322,9122,7043,8901,1,0,,,,,,,,"));
			coordCSV.Add(4324, GetCoordCSV("4324,WGS 72BE,6324,WGS 72 Transit Broadcast Ephemeris,6324,9122,7043,8901,1,0,9606,0,0,1.9,0,0,0.814,-0.38"));
			coordCSV.Add(4326, GetCoordCSV("4326,WGS 84,6326,World Geodetic System 1984,6326,9122,7030,8901,1,0,,,,,,,,"));
			coordCSV.Add(4600, GetCoordCSV("4600,Anguilla 1957,6600,Anguilla 1957,6600,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4601, GetCoordCSV("4601,Antigua 1943,6601,Antigua 1943,6601,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4602, GetCoordCSV("4602,Dominica 1945,6602,Dominica 1945,6602,9122,7012,8901,1,0,9603,725,685,536,,,,"));
			coordCSV.Add(4603, GetCoordCSV("4603,Grenada 1953,6603,Grenada 1953,6603,9122,7012,8901,1,0,9603,72,213.7,93,,,,"));
			coordCSV.Add(4604, GetCoordCSV("4604,Montserrat 1958,6604,Montserrat 1958,6604,9122,7012,8901,1,0,9603,174,359,365,,,,"));
			coordCSV.Add(4605, GetCoordCSV("4605,St. Kitts 1955,6605,St. Kitts 1955,6605,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4606, GetCoordCSV("4606,St. Lucia 1955,6606,St. Lucia 1955,6606,9122,7012,8901,1,0,9603,-149,128,296,,,,"));
			coordCSV.Add(4607, GetCoordCSV("4607,St. Vincent 1945,6607,St. Vincent 1945,6607,9122,7012,8901,1,0,9603,195.671,332.517,274.607,,,,"));
			coordCSV.Add(4608, GetCoordCSV("4608,\"NAD27(76)\",6608,\"North American Datum 1927 (1976)\",6608,9122,7008,8901,1,0,,,,,,,,"));
			coordCSV.Add(4609, GetCoordCSV("4609,\"NAD27(CGQ77)\",6609,\"North American Datum 1927 (CGQ77)\",6609,9122,7008,8901,1,0,,,,,,,,"));
			coordCSV.Add(4610, GetCoordCSV("4610,Xian 1980,6610,Xian 1980,6610,9122,7049,8901,1,0,,,,,,,,"));
			coordCSV.Add(4611, GetCoordCSV("4611,Hong Kong 1980,6611,Hong Kong 1980,6611,9122,7022,8901,1,0,9606,-162.619,-276.959,-161.764,0.067753,-2.243649,-1.158827,-1.094246"));
			coordCSV.Add(4612, GetCoordCSV("4612,JGD2000,6612,Japanese Geodetic Datum 2000,6612,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4613, GetCoordCSV("4613,Segara,6613,Gunung Segara,6613,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4614, GetCoordCSV("4614,QND95,6614,Qatar National Datum 1995,6614,9122,7022,8901,1,0,9606,-119.4248,-303.65872,-11.00061,1.164298,0.174458,1.096259,3.657065"));
			coordCSV.Add(4615, GetCoordCSV("4615,Porto Santo,6615,Porto Santo 1936,6615,9122,7022,8901,1,0,9603,-499,-249,314,,,,"));
			coordCSV.Add(4616, GetCoordCSV("4616,Selvagem Grande,6616,Selvagem Grande,6616,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4617, GetCoordCSV("4617,\"NAD83(CSRS)\",6140,NAD83 Canadian Spatial Reference System,6140,9122,7019,8901,1,0,,,,,,,,"));
			coordCSV.Add(4618, GetCoordCSV("4618,SAD69,6618,South American Datum 1969,6618,9122,7050,8901,1,0,,,,,,,,"));
			coordCSV.Add(4619, GetCoordCSV("4619,SWEREF99,6619,SWEREF99,6619,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4620, GetCoordCSV("4620,Point 58,6620,Point 58,6620,9122,7012,8901,1,0,9603,-106,-129,165,,,,"));
			coordCSV.Add(4621, GetCoordCSV("4621,Fort Marigot,6621,Fort Marigot,6621,9122,7022,8901,1,0,9603,137,248,-430,,,,"));
			coordCSV.Add(4622, GetCoordCSV("4622,Guadeloupe 1948,6622,Guadeloupe 1948,6622,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4623, GetCoordCSV("4623,CSG67,6623,Centre Spatial Guyanais 1967,6623,9122,7022,8901,1,0,9603,-186,230,110,,,,"));
			coordCSV.Add(4624, GetCoordCSV("4624,RGFG95,6624,Reseau Geodesique Francais Guyane 1995,6624,9122,7019,8901,1,0,9603,2,2,-2,,,,"));
			coordCSV.Add(4625, GetCoordCSV("4625,Martinique 1938,6625,Martinique 1938,6625,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4626, GetCoordCSV("4626,Reunion 1947,6626,Reunion 1947,6626,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4627, GetCoordCSV("4627,RGR92,6627,Reseau Geodesique de la Reunion 1992,6627,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4628, GetCoordCSV("4628,Tahiti 52,6628,Tahiti 52,6628,9122,7022,8901,1,0,9603,162,117,154,,,,"));
			coordCSV.Add(4629, GetCoordCSV("4629,Tahaa 54,6629,Tahaa 54,6629,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4630, GetCoordCSV("4630,IGN72 Nuku Hiva,6630,IGN72 Nuku Hiva,6630,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4631, GetCoordCSV("4631,K0 1949,6631,K0 1949,6631,9122,7022,8901,1,1,9603,145,-187,103,,,,"));
			coordCSV.Add(4632, GetCoordCSV("4632,Combani 1950,6632,Combani 1950,6632,9122,7022,8901,1,0,9603,-382,-59,-262,,,,"));
			coordCSV.Add(4633, GetCoordCSV("4633,IGN56 Lifou,6633,IGN56 Lifou,6633,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4634, GetCoordCSV("4634,IGN72 Grand Terre,6634,IGN72 Grande Terre,6634,9108,7022,8901,1,1,,,,,,,,"));
			coordCSV.Add(4635, GetCoordCSV("4635,ST87 Ouvea,6635,ST87 Ouvea,6635,9122,7022,8901,1,1,9606,-122.383,-188.696,103.344,3.5107,-4.9668,-5.7047,4.4798"));
			coordCSV.Add(4636, GetCoordCSV("4636,Petrels 1972,6636,Petrels 1972,6636,9122,7022,8901,1,0,9603,365,194,166,,,,"));
			coordCSV.Add(4637, GetCoordCSV("4637,Perroud 1950,6637,Pointe Geologie Perroud 1950,6637,9122,7022,8901,1,0,9603,325,154,172,,,,"));
			coordCSV.Add(4638, GetCoordCSV("4638,Saint Pierre et Miquelon 1950,6638,Saint Pierre et Miquelon 1950,6638,9122,7008,8901,1,0,9603,30,430,368,,,,"));
			coordCSV.Add(4639, GetCoordCSV("4639,MOP78,6639,MOP78,6639,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4640, GetCoordCSV("4640,RRAF 1991,6640,Reseau de Reference des Antilles Francaises 1991,6640,9122,7030,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4641, GetCoordCSV("4641,IGN53 Mare,6641,IGN53 Mare,6641,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4642, GetCoordCSV("4642,ST84 Ile des Pins,6642,ST84 Ile des Pins,6642,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4643, GetCoordCSV("4643,ST71 Belep,6643,ST71 Belep,6643,9122,7022,8901,1,0,9606,-480.26,-438.32,-643.429,16.3119,20.1721,-4.0349,-111.7002"));
			coordCSV.Add(4644, GetCoordCSV("4644,NEA74 Noumea,6644,NEA74 Noumea,6644,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4645, GetCoordCSV("4645,RGNC 1991,6645,Reseau Geodesique Nouvelle Caledonie 1991,6645,9122,7022,8901,1,1,9603,0,0,0,,,,"));
			coordCSV.Add(4646, GetCoordCSV("4646,Grand Comoros,6646,Grand Comoros,6646,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4657, GetCoordCSV("4657,Reykjavik 1900,6657,Reykjavik 1900,6657,9122,7051,8901,1,0,9603,-28,199,5,,,,"));
			coordCSV.Add(4658, GetCoordCSV("4658,Hjorsey 1955,6658,Hjorsey 1955,6658,9122,7022,8901,1,0,9603,-73,46,-86,,,,"));
			coordCSV.Add(4659, GetCoordCSV("4659,ISN93,6659,Islands Network 1993,6659,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4660, GetCoordCSV("4660,Helle 1954,6660,Helle 1954,6660,9122,7022,8901,1,0,9606,982.6087,552.753,-540.873,32.39344,-153.25684,-96.2266,16.805"));
			coordCSV.Add(4661, GetCoordCSV("4661,LKS92,6661,Latvia 1992,6661,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4662, GetCoordCSV("4662,IGN72 Grande Terre,6634,IGN72 Grande Terre,6634,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4663, GetCoordCSV("4663,Porto Santo 1995,6663,Porto Santo 1995,6663,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4664, GetCoordCSV("4664,Azores Oriental 1995,6664,Azores Oriental Islands 1995,6664,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4665, GetCoordCSV("4665,Azores Central 1995,6665,Azores Central Islands 1995,6665,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4666, GetCoordCSV("4666,Lisbon 1890,6666,Lisbon 1890,6666,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4667, GetCoordCSV("4667,IKBD-92,6667,Iraq-Kuwait Boundary Datum 1992,6667,9122,7030,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4668, GetCoordCSV("4668,ED79,6668,European Datum 1979,6668,9122,7022,8901,1,0,9603,-86,-98,-119,,,,"));
			coordCSV.Add(4669, GetCoordCSV("4669,LKS94,6126,\"Lithuania 1994 (ETRS89)\",6126,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4670, GetCoordCSV("4670,IGM95,6670,Istituto Geografico Militaire 1995,6670,9122,7030,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4671, GetCoordCSV("4671,Voirol 1879,6671,Voirol 1879,6671,9122,7011,8901,1,0,,,,,,,,"));
			coordCSV.Add(4672, GetCoordCSV("4672,CI1971,6672,Chatham Islands Datum 1971,6672,9122,7022,8901,1,0,9603,175,-38,113,,,,"));
			coordCSV.Add(4673, GetCoordCSV("4673,CI1979,6673,Chatham Islands Datum 1979,6673,9122,7022,8901,1,0,9607,174.05,-25.49,112.57,0,0,-0.554,0.2263"));
			coordCSV.Add(4674, GetCoordCSV("4674,SIRGAS 2000,6674,Sistema de Referencia Geocentrico para America del Sur 2000,6674,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4675, GetCoordCSV("4675,Guam 1963,6675,Guam 1963,6675,9122,7008,8901,1,0,9603,-100,-248,259,,,,"));
			coordCSV.Add(4676, GetCoordCSV("4676,Vientiane 1982,6676,Vientiane 1982,6676,9122,7024,8901,1,0,,,,,,,,"));
			coordCSV.Add(4677, GetCoordCSV("4677,Lao 1993,6677,Lao 1993,6677,9122,7024,8901,1,0,,,,,,,,"));
			coordCSV.Add(4678, GetCoordCSV("4678,Lao 1997,6678,Lao National Datum 1997,6678,9122,7024,8901,1,0,9603,44.585,-131.212,-39.544,,,,"));
			coordCSV.Add(4679, GetCoordCSV("4679,Jouik 1961,6679,Jouik 1961,6679,9122,7012,8901,1,0,9603,-80.01,253.26,291.19,,,,"));
			coordCSV.Add(4680, GetCoordCSV("4680,Nouakchott 1965,6680,Nouakchott 1965,6680,9122,7012,8901,1,0,9603,124.5,-63.5,-281,,,,"));
			coordCSV.Add(4681, GetCoordCSV("4681,Mauritania 1999,6681,Mauritania 1999,6681,9122,7012,8901,1,1,,,,,,,,"));
			coordCSV.Add(4682, GetCoordCSV("4682,Gulshan 303,6682,Gulshan 303,6682,9122,7015,8901,1,0,,,,,,,,"));
			coordCSV.Add(4683, GetCoordCSV("4683,PRS92,6683,Philippine Reference System 1992,6683,9122,7008,8901,1,0,9607,-127.62,-67.24,-47.04,3.068,-4.903,-1.578,-1.06"));
			coordCSV.Add(4684, GetCoordCSV("4684,Gan 1970,6684,Gan 1970,6684,9122,7022,8901,1,0,9603,-133,-321,50,,,,"));
			coordCSV.Add(4685, GetCoordCSV("4685,Gandajika,6685,Gandajika,6685,9122,7022,8901,1,1,,,,,,,,"));
			coordCSV.Add(4686, GetCoordCSV("4686,MAGNA-SIRGAS,6686,Marco Geocentrico Nacional de Referencia,6686,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4687, GetCoordCSV("4687,RGPF,6687,Reseau Geodesique de la Polynesie Francaise,6687,9122,7019,8901,1,0,,,,,,,,"));
			coordCSV.Add(4688, GetCoordCSV("4688,Fatu Iva 72,6688,Fatu Iva 72,6688,9122,7022,8901,1,0,9607,347.103,1078.125,2623.922,33.8875,-70.6773,9.3943,186.074"));
			coordCSV.Add(4689, GetCoordCSV("4689,IGN63 Hiva Oa,6689,IGN63 Hiva Oa,6689,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4690, GetCoordCSV("4690,Tahiti 79,6690,Tahiti 79,6690,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4691, GetCoordCSV("4691,Moorea 87,6691,Moorea 87,6691,9122,7022,8901,1,0,9607,215.525,149.593,176.229,3.2624,1.692,1.1571,10.4773"));
			coordCSV.Add(4692, GetCoordCSV("4692,Maupiti 83,6692,Maupiti 83,6692,9122,7022,8901,1,0,9603,217.037,86.959,23.956,,,,"));
			coordCSV.Add(4693, GetCoordCSV("4693,Nakhl-e Ghanem,6693,Nakhl-e Ghanem,6693,9122,7030,8901,1,0,9603,0,-0.15,0.68,,,,"));
			coordCSV.Add(4694, GetCoordCSV("4694,POSGAR 94,6694,Posiciones Geodesicas Argentinas 1994,6694,9122,7019,8901,1,0,,,,,,,,"));
			coordCSV.Add(4695, GetCoordCSV("4695,Katanga 1955,6695,Katanga 1955,6695,9122,7008,8901,1,0,,,,,,,,"));
			coordCSV.Add(4696, GetCoordCSV("4696,Kasai 1953,6696,Kasai 1953,6696,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4697, GetCoordCSV("4697,IGC 1962 6th Parallel South,6697,IGC 1962 Arc of the 6th Parallel South,6697,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4698, GetCoordCSV("4698,IGN 1962 Kerguelen,6698,IGN 1962 Kerguelen,6698,9122,7022,8901,1,0,9603,145,-187,103,,,,"));
			coordCSV.Add(4699, GetCoordCSV("4699,Le Pouce 1934,6699,Le Pouce 1934,6699,9122,7012,8901,1,0,9603,-770.1,158.4,-498.2,,,,"));
			coordCSV.Add(4700, GetCoordCSV("4700,IGN Astro 1960,6700,IGN Astro 1960,6700,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4701, GetCoordCSV("4701,IGCB 1955,6701,Institut Geographique du Congo Belge 1955,6701,9122,7012,8901,1,0,9603,-79.9,-158,-168.9,,,,"));
			coordCSV.Add(4702, GetCoordCSV("4702,Mauritania 1999,6702,Mauritania 1999,6702,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4703, GetCoordCSV("4703,Mhast 1951,6703,Missao Hidrografico Angola y Sao Tome 1951,6703,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4704, GetCoordCSV("4704,\"Mhast (onshore)\",6704,\"Mhast (onshore)\",6704,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4705, GetCoordCSV("4705,\"Mhast (offshore)\",6705,\"Mhast (offshore)\",6705,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4706, GetCoordCSV("4706,Egypt Gulf of Suez S-650 TL,6706,Egypt Gulf of Suez S-650 TL,6706,9122,7020,8901,1,0,9603,-146.21,112.63,4.05,,,,"));
			coordCSV.Add(4707, GetCoordCSV("4707,Tern Island 1961,6707,Tern Island 1961,6707,9122,7022,8901,1,0,9603,114,-116,-333,,,,"));
			coordCSV.Add(4708, GetCoordCSV("4708,Cocos Islands 1965,6708,Cocos Islands 1965,6708,9122,7003,8901,1,0,9603,-491,-22,435,,,,"));
			coordCSV.Add(4709, GetCoordCSV("4709,Iwo Jima 1945,6709,Iwo Jima 1945,6709,9122,7022,8901,1,0,9603,145,75,-272,,,,"));
			coordCSV.Add(4710, GetCoordCSV("4710,St. Helena 1971,6710,St. Helena 1971,6710,9122,7022,8901,1,0,9603,-320,550,-494,,,,"));
			coordCSV.Add(4711, GetCoordCSV("4711,Marcus Island 1952,6711,Marcus Island 1952,6711,9122,7022,8901,1,0,9603,124,-234,-25,,,,"));
			coordCSV.Add(4712, GetCoordCSV("4712,Ascension Island 1958,6712,Ascension Island 1958,6712,9122,7022,8901,1,0,9603,-205,107,53,,,,"));
			coordCSV.Add(4713, GetCoordCSV("4713,Ayabelle Lighthouse,6713,Ayabelle Lighthouse,6713,9122,7012,8901,1,0,9603,-79,-129,145,,,,"));
			coordCSV.Add(4714, GetCoordCSV("4714,Bellevue,6714,Bellevue,6714,9122,7022,8901,1,0,9603,-127,-769,472,,,,"));
			coordCSV.Add(4715, GetCoordCSV("4715,Camp Area Astro,6715,Camp Area Astro,6715,9122,7022,8901,1,0,9603,-104,-129,239,,,,"));
			coordCSV.Add(4716, GetCoordCSV("4716,Phoenix Islands 1966,6716,Phoenix Islands 1966,6716,9122,7022,8901,1,0,9603,298,-304,-375,,,,"));
			coordCSV.Add(4717, GetCoordCSV("4717,Cape Canaveral,6717,Cape Canaveral,6717,9122,7008,8901,1,0,9603,-2,151,181,,,,"));
			coordCSV.Add(4718, GetCoordCSV("4718,Solomon 1968,6718,Solomon 1968,6718,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4719, GetCoordCSV("4719,Easter Island 1967,6719,Easter Island 1967,6719,9122,7022,8901,1,0,9603,211,147,111,,,,"));
			coordCSV.Add(4720, GetCoordCSV("4720,Fiji 1986,6720,Fiji Geodetic Datum 1986,6720,9122,7043,8901,1,0,,,,,,,,"));
			coordCSV.Add(4721, GetCoordCSV("4721,Fiji 1956,6721,Fiji 1956,6721,9122,7022,8901,1,0,9603,265.025,384.929,-194.046,,,,"));
			coordCSV.Add(4722, GetCoordCSV("4722,South Georgia 1968,6722,South Georgia 1968,6722,9122,7022,8901,1,0,9603,-794,119,-298,,,,"));
			coordCSV.Add(4723, GetCoordCSV("4723,Grand Cayman 1959,6723,Grand Cayman 1959,6723,9122,7008,8901,1,0,9603,67.8,106.1,138.8,,,,"));
			coordCSV.Add(4724, GetCoordCSV("4724,Diego Garcia 1969,6724,Diego Garcia 1969,6724,9122,7022,8901,1,0,9603,208,-435,-229,,,,"));
			coordCSV.Add(4725, GetCoordCSV("4725,Johnston Island 1961,6725,Johnston Island 1961,6725,9122,7022,8901,1,0,9603,189,-79,-202,,,,"));
			coordCSV.Add(4726, GetCoordCSV("4726,Little Cayman 1961,6726,Little Cayman 1961,6726,9122,7008,8901,1,0,,,,,,,,"));
			coordCSV.Add(4727, GetCoordCSV("4727,Midway 1961,6727,Midway 1961,6727,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4728, GetCoordCSV("4728,Pico de la Nieves,6728,Pico de la Nieves,6728,9122,7022,8901,1,0,9603,-307,-92,127,,,,"));
			coordCSV.Add(4729, GetCoordCSV("4729,Pitcairn 1967,6729,Pitcairn 1967,6729,9122,7022,8901,1,0,9603,185,165,42,,,,"));
			coordCSV.Add(4730, GetCoordCSV("4730,Santo 1965,6730,Santo 1965,6730,9122,7022,8901,1,0,9603,170,42,84,,,,"));
			coordCSV.Add(4731, GetCoordCSV("4731,Viti Levu 1916,6731,Viti Levu 1916,6731,9122,7012,8901,1,1,9603,51,391,-36,,,,"));
			coordCSV.Add(4732, GetCoordCSV("4732,Marshall Islands 1960,6732,Marshall Islands 1960,6732,9122,7053,8901,1,0,9603,102,52,-38,,,,"));
			coordCSV.Add(4733, GetCoordCSV("4733,Wake Island 1952,6733,Wake Island 1952,6733,9122,7022,8901,1,0,9603,276,-57,149,,,,"));
			coordCSV.Add(4734, GetCoordCSV("4734,Tristan 1968,6734,Tristan 1968,6734,9122,7022,8901,1,0,9603,-632,438,-609,,,,"));
			coordCSV.Add(4735, GetCoordCSV("4735,Kusaie 1951,6735,Kusaie 1951,6735,9122,7022,8901,1,0,9603,647,1777,-1124,,,,"));
			coordCSV.Add(4736, GetCoordCSV("4736,Deception Island,6736,Deception Island,6736,9122,7012,8901,1,0,9603,260,12,-147,,,,"));
			coordCSV.Add(4737, GetCoordCSV("4737,Korea 2000,6737,Geocentric datum of Korea,6737,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4738, GetCoordCSV("4738,Hong Kong 1963,6738,Hong Kong 1963,6738,9122,7007,8901,1,0,,,,,,,,"));
			coordCSV.Add(4739, GetCoordCSV("4739,\"Hong Kong 1963(67)\",6739,\"Hong Kong 1963(67)\",6739,9122,7022,8901,1,0,9606,-156,-271,-189,,,,"));
			coordCSV.Add(4740, GetCoordCSV("4740,PZ-90,6740,Parametrop Zemp 1990,6740,9122,7054,8901,1,0,9607,0,0,1.5,0,0,-0.076,0"));
			coordCSV.Add(4741, GetCoordCSV("4741,FD54,6741,Faroe Datum 1954,6741,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4742, GetCoordCSV("4742,GDM2000,6742,Geodetic Datum of Malaysia 2000,6742,9122,7019,8901,1,0,,,,,,,,"));
			coordCSV.Add(4743, GetCoordCSV("4743,\"Karbala 1979 (Polservice)\",6743,\"Karbala 1979 (Polservice)\",6743,9122,7012,8901,1,0,9603,84.1,-320.1,218.7,,,,"));
			coordCSV.Add(4744, GetCoordCSV("4744,Nahrwan 1934,6744,Nahrwan 1934,6744,9122,7012,8901,1,0,,,,,,,,"));
			coordCSV.Add(4745, GetCoordCSV("4745,\"RD/83\",6745,\"Rauenberg Datum/83\",6745,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4746, GetCoordCSV("4746,\"PD/83\",6746,\"Potsdam Datum/83\",6746,9122,7004,8901,1,0,,,,,,,,"));
			coordCSV.Add(4747, GetCoordCSV("4747,GR96,6747,Greenland 1996,6747,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4748, GetCoordCSV("4748,Vanua Levu 1915,6748,Vanua Levu 1915,6748,9122,7055,8901,1,0,9603,51,391,-36,,,,"));
			coordCSV.Add(4749, GetCoordCSV("4749,RGNC91-93,6749,Reseau Geodesique de Nouvelle Caledonie 91-93,6749,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4750, GetCoordCSV("4750,ST87 Ouvea,6750,ST87 Ouvea,6750,9122,7030,8901,1,0,9603,-56.263,16.136,-22.856,,,,"));
			coordCSV.Add(4751, GetCoordCSV("4751,\"Kertau (RSO)\",6751,\"Kertau (RSO)\",6751,9122,7056,8901,1,0,,,,,,,,"));
			coordCSV.Add(4752, GetCoordCSV("4752,Viti Levu 1912,6752,Viti Levu 1912,6752,9122,7055,8901,1,0,9603,51,391,-36,,,,"));
			coordCSV.Add(4753, GetCoordCSV("4753,fk89,6753,fk89,6753,9122,7022,8901,1,0,,,,,,,,"));
			coordCSV.Add(4754, GetCoordCSV("4754,LGD2006,6754,Libyan Geodetic Datum 2006,6754,9122,7022,8901,1,0,9603,-208.4058,-109.8777,-2.5764,,,,"));
			coordCSV.Add(4755, GetCoordCSV("4755,DGN95,6755,Datum Geodesi Nasional 1995,6755,9122,7030,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4756, GetCoordCSV("4756,VN-2000,6756,Vietnam 2000,6756,9122,7030,8901,1,0,,,,,,,,"));
			coordCSV.Add(4757, GetCoordCSV("4757,SVY21,6757,SVY21,6757,9122,7030,8901,1,0,,,,,,,,"));
			coordCSV.Add(4758, GetCoordCSV("4758,JAD2001,6758,Jamaica 2001,6758,9122,7030,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4759, GetCoordCSV("4759,\"NAD83(NSRS2007)\",6759,\"NAD83 (National Spatial Reference System 2007)\",6759,9122,7019,8901,1,0,9603,0,0,0,,,,"));
			coordCSV.Add(4760, GetCoordCSV("4760,WGS 66,6760,World Geodetic System 1966,6760,9122,7025,8901,1,0,,,,,,,,"));
			coordCSV.Add(4801, GetCoordCSV("4801,\"Bern 1898 (Bern)\",6801,\"CH1903 (Bern)\",6149,9122,7004,8907,1,0,,,,,,,,"));
			coordCSV.Add(4802, GetCoordCSV("4802,\"Bogota 1975 (Bogota)\",6802,\"Bogota 1975 (Bogota)\",6218,9122,7022,8904,1,0,,,,,,,,"));
			coordCSV.Add(4803, GetCoordCSV("4803,\"Lisbon (Lisbon)\",6803,\"Lisbon 1937 (Lisbon)\",6207,9122,7022,8902,1,0,,,,,,,,"));
			coordCSV.Add(4804, GetCoordCSV("4804,\"Makassar (Jakarta)\",6804,\"Makassar (Jakarta)\",6257,9122,7004,8908,1,0,9603,-587.8,519.75,145.76,,,,"));
			coordCSV.Add(4805, GetCoordCSV("4805,\"MGI (Ferro)\",6805,\"Militar-Geographische Institut (Ferro)\",6312,9122,7004,8909,1,0,,,,,,,,"));
			coordCSV.Add(4806, GetCoordCSV("4806,\"Monte Mario (Rome)\",6806,\"Monte Mario (Rome)\",6265,9122,7022,8906,1,0,,,,,,,,"));
			coordCSV.Add(4807, GetCoordCSV("4807,\"NTF (Paris)\",6807,\"Nouvelle Triangulation Francaise (Paris)\",6275,9105,7011,8903,1,0,9603,-168,-60,320,,,,"));
			coordCSV.Add(4808, GetCoordCSV("4808,\"Padang (Jakarta)\",6808,\"Padang 1884 (Jakarta)\",6280,9122,7004,8908,1,0,,,,,,,,"));
			coordCSV.Add(4809, GetCoordCSV("4809,\"Belge 1950 (Brussels)\",6809,\"Reseau National Belge 1950 (Brussels)\",6215,9122,7022,8910,1,0,,,,,,,,"));
			coordCSV.Add(4810, GetCoordCSV("4810,\"Tananarive (Paris)\",6810,\"Tananarive 1925 (Paris)\",6297,9105,7022,8903,1,0,9603,-189,-242,-91,,,,"));
			coordCSV.Add(4811, GetCoordCSV("4811,\"Voirol 1875 (Paris)\",6811,\"Voirol 1875 (Paris)\",6304,9105,7011,8903,1,0,9603,-73,-247,227,,,,"));
			coordCSV.Add(4813, GetCoordCSV("4813,\"Batavia (Jakarta)\",6813,\"Batavia (Jakarta)\",6211,9122,7004,8908,1,0,,,,,,,,"));
			coordCSV.Add(4814, GetCoordCSV("4814,\"RT38 (Stockholm)\",6814,\"Stockholm 1938 (Stockholm)\",6308,9122,7004,8911,1,0,,,,,,,,"));
			coordCSV.Add(4815, GetCoordCSV("4815,\"Greek (Athens)\",6815,\"Greek (Athens)\",6120,9122,7004,8912,1,0,,,,,,,,"));
			coordCSV.Add(4816, GetCoordCSV("4816,\"Carthage (Paris)\",6816,\"Carthage (Paris)\",6223,9105,7011,8903,1,0,,,,,,,,"));
			coordCSV.Add(4817, GetCoordCSV("4817,\"NGO 1948 (Oslo)\",6817,\"NGO 1948 (Oslo)\",6273,9122,7005,8913,1,0,9606,278.3,93,474.5,7.889,0.05,-6.61,6.21"));
			coordCSV.Add(4818, GetCoordCSV("4818,\"S-JTSK (Ferro)\",6818,\"S-JTSK (Ferro)\",6156,9122,7004,8909,1,0,,,,,,,,"));
			coordCSV.Add(4819, GetCoordCSV("4819,\"Nord Sahara 1959 (Paris)\",6819,\"Nord Sahara 1959 (Paris)\",6307,9105,7012,8903,1,1,,,,,,,,"));
			coordCSV.Add(4820, GetCoordCSV("4820,\"Segara (Jakarta)\",6820,\"Gunung Segara (Jakarta)\",6613,9122,7004,8908,1,0,,,,,,,,"));
			coordCSV.Add(4821, GetCoordCSV("4821,\"Voirol 1879 (Paris)\",6821,\"Voirol 1879 (Paris)\",6821,9105,7011,8903,1,0,,,,,,,,"));
			coordCSV.Add(4901, GetCoordCSV("4901,\"ATF (Paris)\",6901,\"Ancienne Triangulation Francaise (Paris)\",6901,9105,7027,8903,1,0,,,,,,,,"));
			coordCSV.Add(4902, GetCoordCSV("4902,\"NDG (Paris)\",6902,\"Nord de Guerre (Paris)\",6902,9105,7027,8903,1,0,,,,,,,,"));
			coordCSV.Add(4903, GetCoordCSV("4903,\"Madrid 1870 (Madrid)\",6903,\"Madrid 1870 (Madrid)\",6903,9122,7028,8905,1,0,,,,,,,,"));
			coordCSV.Add(4904, GetCoordCSV("4904,\"Lisbon 1890 (Lisbon)\",6904,\"Lisbon 1890 (Lisbon)\",6666,9122,7004,8902,1,0,,,,,,,,"));
			#endregion
		}
	}
}
