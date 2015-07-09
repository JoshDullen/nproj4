//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Initialize projection object from string definition. Includes
//			pj_init(), pj_init_plus() and pj_free() function.
// Author:	Gerald Evenden, Frank Warmerdam <warmerdam@pobox.com>
//
//*****************************************************************************
// Copyright (c) 1995, Gerald Evenden
// Copyright (c) 2002, Frank Warmerdam <warmerdam@pobox.com>
// Copyright (c) 2009-2011 by the Authors
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
using System.IO;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		//**********************************************************************
		//								get_opt()
		//**********************************************************************
		static string get_opt(projCtx ctx, List<string> start, StreamReader fid, string name, out bool found_def)
		{
			bool first=true, done=false;
			string next=null;

			found_def=false;

			while(!fid.EndOfStream&&!done)
			{
				string line=fid.ReadLine();
				if(line.Length<2) continue;
				line=line.TrimStart();
				if(line.Length<2) continue;
				if(line[0]=='#') continue;

				string[] words=line.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
				if(words.Length<1) continue;

				foreach(string word in words)
				{
					if(word[0]=='#') break;
					if(word[0]=='<')
					{ // control name
						if(!first) { done=true; break; }
						if(first&&word=="<"+name+">")
						{
							first=false;
							found_def=true;
						}
						continue;
					}
					else if(!first&&!pj_param_t(ctx, start, word))
					{
						// don't default ellipse if datum, ellps or any earth model information is set.
						if(word.StartsWith("ellps=")||
							(!pj_param_t(ctx, start, "datum")&&!pj_param_t(ctx, start, "ellps")&&!pj_param_t(ctx, start, "a")&&
							!pj_param_t(ctx, start, "b")&&!pj_param_t(ctx, start, "rf")&&!pj_param_t(ctx, start, "f")))
						{
							next=pj_mkparam(word);
							if(next==null) return null;
							start.Add(next);
						}
					}
				}
			}

			return next;
		}

		//**********************************************************************
		//							get_defaults()
		//**********************************************************************
		static string get_defaults(projCtx ctx, List<string> start, string name)
		{
			string next=null;
			try
			{
				FileStream fs=pj_open_lib(ctx, "proj_def.dat", FileAccess.Read);
				if(fs!=null)
				{
					StreamReader fid=new StreamReader(fs);
					bool found_def;
					next=get_opt(ctx, start, fid, "general", out found_def);
					fid.BaseStream.Seek(0, SeekOrigin.Begin);
					next=get_opt(ctx, start, fid, name, out found_def);
					fid.Close();
				}
				else
				{
					// if proj_def.dat not found
					start.Add(pj_mkparam("ellps=WGS84"));

					if(name=="aea")
					{
						start.Add(pj_mkparam("lat_1=29.5"));
						start.Add(pj_mkparam("lat_2=45.5"));
					}
					else if(name=="lcc")
					{
						start.Add(pj_mkparam("lat_1=33"));
						start.Add(pj_mkparam("lat_2=45"));
					}
					else if(name=="lagrng") start.Add(pj_mkparam("W=2"));
				}
			}
			catch
			{
			}

			if(Libc.errno!=0) Libc.errno=0; // don't care if can't open file
			ctx.last_errno=0;

			return next;
		}

		//**********************************************************************
		//								get_init()
		//**********************************************************************
		static string get_init(projCtx ctx, List<string> start, string name, out bool found_def)
		{
			found_def=false;

			int opt=name.LastIndexOf(':');
			if(opt==-1) { pj_ctx_set_errno(ctx, -3); return null; }
			if(opt==name.Length-1) { pj_ctx_set_errno(ctx, -3); return null; } // last char ':' => no name given

			string next=null;
			try
			{
				FileStream fs=pj_open_lib(ctx, name.Substring(0, opt), FileAccess.Read);
				if(fs==null) return null;

				StreamReader fid=new StreamReader(fs);
				next=get_opt(ctx, start, fid, name.Substring(opt+1), out found_def);
				fid.Close();
			}
			catch
			{
			}

			return next;
		}

		//**********************************************************************
		//							pj_init_plus()
		//
		//		Same as pj_init() except it takes one argument string with
		//		individual arguments preceeded by '+', such as "+proj=utm
		//		+zone=11 +ellps=WGS84".
		//**********************************************************************
		public static PJ pj_init_plus(string definition)
		{
			return pj_init_plus_ctx(pj_get_default_ctx(), definition);
		}

		public static PJ pj_init_plus_ctx(projCtx ctx, string definition)
		{
			definition=definition.Trim(' ', '\t', '\n', '\r');
			if(definition.Length<=0) { pj_ctx_set_errno(ctx, -1); return null; }

			// split on whitespace not in double-quotes-strings (where \" is a double-quote in a string)
			List<string> argv=new List<string>();

			StringBuilder part=new StringBuilder();
			bool inString=false;
			bool wasBackslash=false;
			for(int i=0; i<definition.Length; i++)
			{
				if(wasBackslash)
				{
					wasBackslash=false;

					if(definition[i]=='"')
					{
						part.Append('"');
						continue;
					}

					part.Append('\\');
				}

				switch(definition[i])
				{
					case '"': inString=!inString; break;
					case '\\': wasBackslash=true; break;
					case ' ':
					case '\t':
					case '\n':
					case '\r':
						if(inString) part.Append(definition[i]);
						else
						{
							string str=part.ToString().Trim();
							if(str.Length>1&&str.StartsWith("+")) argv.Add(str.Substring(1));

							part.Clear();
						}
						break;
					default: part.Append(definition[i]); break;
				}
			}

			string strLast=part.ToString().Trim();
			if(strLast.Length>1&&strLast.StartsWith("+")) argv.Add(strLast.Substring(1));

			// perform actual initialization
			return pj_init_ctx(ctx, argv.ToArray());
		}

		//**********************************************************************
		//								pj_init()
		//
		//		Main entry point for initialing a PJ projections
		//		definition. Note that the projection specific function is
		//		called to do the initial allocation so it can be created
		//		large enough to hold projection specific parameters.
		//**********************************************************************
		public static PJ pj_init(string[] argv)
		{
			return pj_init_ctx(pj_get_default_ctx(), argv);
		}

		public static PJ pj_init_ctx(projCtx ctx, string[] argv)
		{
			ctx.last_errno=0;

			if(argv.Length<=0) { pj_ctx_set_errno(ctx, -1); return null; }

			List<string> start=new List<string>();

			Libc.errno=pj_errno=0;
			
			// put arguments into internal linked list
			for(int a=0; a<argv.Length; a++)
			{
				string curr=pj_mkparam(argv[a]);
				if(curr==null) { pj_ctx_set_errno(ctx, Libc.errno); return null; }
				start.Add(curr);
			}

			// check if +init present
			if(pj_param_t(ctx, start, "init"))
			{
				bool found_def;
				string curr=get_init(ctx, start, pj_param_s(ctx, start, "init"), out found_def);
				if(curr==null||!found_def)
				{
					if(pj_errno!=0||Libc.errno!=0)
					{
						if(pj_errno==0) pj_ctx_set_errno(ctx, Libc.errno);
					}
					else pj_ctx_set_errno(ctx, -2);

					return null;
				}
			}

			// find projection selection
			string name=pj_param_s(ctx, start, "proj");
			if(name=="") { pj_ctx_set_errno(ctx, -4); return null; }

			// set defaults, unless inhibited
			if(!pj_param_b(ctx, start, "no_defs")) get_defaults(ctx, start, name);

			// allocate projection structure
			PJ PIN=GetPJ(name);
			if(PIN==null) { pj_ctx_set_errno(ctx, -5); return null; }
			PIN.ctx=ctx;
			PIN.parameters=start;
			PIN.is_latlong=false;
			PIN.is_geocent=false;
			PIN.is_long_wrap_set=false;
			PIN.long_wrap_center=0.0;
			PIN.axis="enu";

			PIN.gridlist=null;
			PIN.vgridlist_geoid=null;

			// set datum parameters
			if(pj_datum_set(ctx, start, PIN)) return null;

			// set ellipsoid/sphere parameters
			if(pj_ell_set(ctx, start, out PIN.a, out PIN.es)) return null;

			PIN.a_orig=PIN.a;
			PIN.es_orig=PIN.es;

			PIN.e=Math.Sqrt(PIN.es);
			PIN.ra=1.0/PIN.a;
			PIN.one_es=1.0-PIN.es;
			if(PIN.one_es==0.0) { pj_ctx_set_errno(ctx, -6); return null; }
			PIN.rone_es=1.0/PIN.one_es;

			// Now that we have ellipse information check for WGS84 datum
			if(PIN.datum_type==PJD._3PARAM&&PIN.datum_params[0]==0.0&&PIN.datum_params[1]==0.0&&PIN.datum_params[2]==0.0&&
				PIN.a==6378137.0&&Math.Abs(PIN.es-0.006694379990)<0.000000000050) PIN.datum_type=PJD.WGS84; //WGS84/GRS80

			// set PIN.geoc coordinate system
			PIN.geoc=PIN.es!=0.0&&pj_param_b(ctx, start, "geoc");

			// over-ranging flag
			PIN.over=pj_param_b(ctx, start, "over");

			// vertical datum geoid grids
			PIN.has_geoid_vgrids=pj_param_t(ctx, start, "geoidgrids");
			if(PIN.has_geoid_vgrids) // we need to mark it as used.
				pj_param_s(ctx, start, "geoidgrids");

			// longitude center for wrapping
			PIN.is_long_wrap_set=pj_param_b(ctx, start, "lon_wrap");
			if(PIN.is_long_wrap_set)
				PIN.long_wrap_center=pj_param_r(ctx, start, "lon_wrap");

			// axis orientation
			string axis_arg=pj_param_s(ctx, start,"axis");
			if(axis_arg!=null&&axis_arg.Length!=0)
			{
				string axis_legal="ewnsud";
				if(axis_arg.Length!=3)
				{
					pj_ctx_set_errno(ctx, (int)PJD_ERR.AXIS);
					return null;
				}

				if(axis_legal.IndexOf(axis_arg[0])==-1||axis_legal.IndexOf(axis_arg[1])==-1||(axis_arg.Length>=3&&axis_legal.IndexOf(axis_arg[2])==-1))
				{
					pj_ctx_set_errno(ctx, (int)PJD_ERR.AXIS);
					return null;
				}

				// it would be nice to validate we don't have on axis repeated
				PIN.axis=axis_arg;
			}

			// central meridian
			PIN.lam0=pj_param_r(ctx, start, "lon_0");

			// central latitude
			PIN.phi0=pj_param_r(ctx, start, "lat_0");

			// false easting and northing
			PIN.x0=pj_param_d(ctx, start, "x_0");
			PIN.y0=pj_param_d(ctx, start, "y_0");

			// general scaling factor
			if(pj_param_t(ctx, start, "k_0")) PIN.k0=pj_param_d(ctx, start, "k_0");
			else if(pj_param_t(ctx, start, "k")) PIN.k0=pj_param_d(ctx, start, "k");
			else PIN.k0=1.0;

			if(PIN.k0<=0.0) { pj_ctx_set_errno(ctx, -31); return null; }

			// set units
			double to_meter=double.NaN;
			name=pj_param_s(ctx, start, "units");
			if(name!="")
			{
				to_meter=GetUnitFactor(name);
				if(double.IsNaN(to_meter)) { pj_ctx_set_errno(ctx, -7); return null; }
				PIN.to_meter=to_meter;
			}

			if(double.IsNaN(to_meter))
			{
				string s=pj_param_s(ctx, start, "to_meter");
				if(s!="")
				{
					PIN.to_meter=Libc.strtod(s, out s);
					if(s.Length>0&&s[0]=='/') PIN.to_meter/=Libc.strtod(s.Substring(1), out s); // ratio number
					PIN.fr_meter=1.0/PIN.to_meter;
				}
				else PIN.to_meter=PIN.fr_meter=1.0;
			}
			else PIN.fr_meter=1.0/PIN.to_meter;

			// set vertical units
			double vto_meter=double.NaN;
			name=pj_param_s(ctx, start, "vunits");
			if(name!="")
			{
				vto_meter=GetUnitFactor(name);
				if(double.IsNaN(vto_meter)) { pj_ctx_set_errno(ctx, -7); return null; }
				PIN.vto_meter=vto_meter;
			}

			if(double.IsNaN(vto_meter))
			{
				string s=pj_param_s(ctx, start, "vto_meter");
				if(s!="")
				{
					PIN.vto_meter=Libc.strtod(s, out s);
					if(s.Length>0&&s[0]=='/') PIN.vto_meter/=Libc.strtod(s.Substring(1), out s); // ratio number
					PIN.vfr_meter=1.0/PIN.vto_meter;
				}
				else
				{
					PIN.vto_meter=PIN.to_meter;
					PIN.vfr_meter=PIN.fr_meter;
				}
			}
			else PIN.vfr_meter=1.0/PIN.vto_meter;

			// prime meridian
			name=pj_param_s(ctx, start, "pm");
			if(name!="")
			{
				string value=GetPrimeMeridian(name);

				string next_str;
				if(value==null)
				{
					double tmp=dmstor_ctx(ctx, name, out next_str);
					if((tmp!=0.0||name[0]=='0')&&next_str=="") value=name;
				}

				if(value==null) { pj_ctx_set_errno(ctx, -46); return null; }
				PIN.from_greenwich=dmstor_ctx(ctx, value, out next_str);
			}
			else PIN.from_greenwich=0.0;

			// projection specific initialization
			PIN=PIN.Init();
			if(PIN==null||ctx.last_errno!=0)
			{
				// cleanup error return
				return null;
			}

			return PIN;
		}
	}
}
