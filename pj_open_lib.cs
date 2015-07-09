//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Implementation of pj_open_lib(), and pj_set_finder(). These
//			provide a standard interface for opening projections support
//			data files.
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
using System.IO;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		public delegate string pj_finder_proc(string id); 
		static pj_finder_proc pj_finder=null;
		static string[] search_path=null;

		//**********************************************************************
		//							pj_set_finder()
		//**********************************************************************
		public static void pj_set_finder(pj_finder_proc new_finder)
		{
			pj_finder=new_finder;
		}

		//**********************************************************************
		//							pj_set_searchpath()
		//
		//		Path control for callers that can't practically provide
		//		pj_set_finder() style callbacks. Call with (null) as args
		//		to clear the searchpath set.
		//**********************************************************************
		public static void pj_set_searchpath(string[] path)
		{
			search_path=null;

			if(path==null||path.Length==0) return;

			search_path=new string[path.Length];
			for(int i=0; i<path.Length; i++) search_path[i]=path[i];
		}

		//**********************************************************************
		//							pj_open_lib()
		//**********************************************************************
		public static FileStream pj_open_lib(projCtx ctx, string name, FileAccess mode)
		{
			string sysname;
			string dir_chars="/\\";

			if(name==null||name.Length<3) return null;

			// check if ~/name
			if(name[0]=='~'&&dir_chars.IndexOf(name[1])!=-1)
			{
				sysname=Environment.GetEnvironmentVariable("HOME");
				if(sysname!=null) sysname+=DIR_CHAR+name.Substring(1);
				else return null;
			}
			else
			{
				// or fixed path: "/name", "./name", "../name" or "?:\name"
				if(dir_chars.IndexOf(name[0])!=-1||(name[0]=='.'&&dir_chars.IndexOf(name[1])!=-1)||
					(name[0]=='.'&&name[1]=='.'&&dir_chars.IndexOf(name[2])!=-1)||
					(name[1]==':'&&dir_chars.IndexOf(name[2])!=-1)) sysname=name;
				else
				{
					// or try to use application provided file finder
					if(pj_finder!=null&&pj_finder(name)!=null) sysname=pj_finder(name);
					else
					{
						// or is environment PROJ_LIB defined
						sysname=Environment.GetEnvironmentVariable(PROJ_LIB);
						if(sysname!=null) sysname+=DIR_CHAR+name;
						else
						{
							try
							{
								FileInfo fi=new FileInfo(System.Reflection.Assembly.GetEntryAssembly().Location);
								if(fi.DirectoryName.Length==3) sysname=fi.DirectoryName+PROJ_LIB+DIR_CHAR+name;
								else sysname=fi.DirectoryName+DIR_CHAR+PROJ_LIB+DIR_CHAR+name;
							}
							catch
							{
								sysname=PROJ_LIB+DIR_CHAR+name;
								//else sysname=name; // just try it bare bones
							}
						}
					}
				}
			}

			FileStream fid=null;
			try
			{
				fid=File.Open(sysname, FileMode.Open, mode);
			}
			catch
			{
				fid=null;
			}
			if(fid!=null) Libc.errno=0;

			// If none of those work and we have a search path, try it
			if(fid==null&&search_path!=null)
			{
				foreach(string path in search_path)
				{
					sysname=path+DIR_CHAR+name;
					try
					{
						fid=File.Open(sysname, FileMode.Open, mode);
					}
					catch
					{
						fid=null;
					}
					if(fid!=null) break;
				}
				if(fid!=null) Libc.errno=0;
			}

			if(ctx.last_errno==0&&Libc.errno!=0)
				pj_ctx_set_errno(ctx, Libc.errno);

#if DEBUG
				Console.Error.WriteLine("pj_open_lib({0}): call fopen({1}) - {2}", name, sysname, fid==null?"failed":"succeeded");
#endif

			return fid;
		}
	}
}
