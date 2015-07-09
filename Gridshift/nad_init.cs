//*****************************************************************************
//
// Project: PROJ.4
// Purpose: Load datum shift files into memory.
// Author: Frank Warmerdam, warmerdam@pobox.com
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

using System;
using System.IO;
using System.Text;

namespace Free.Ports.Proj4.Gridshift
{
	public static partial class Grid
	{
		//**********************************************************************
		//							nad_ctable_load()
		//
		//		Load the data portion of a ctable formatted grid.
		//**********************************************************************
		public static bool nad_ctable_load(projCtx ctx, CTABLE ct, Stream fid)
		{
			try
			{
				fid.Seek(80+16+16+8+4, SeekOrigin.Begin); // 80+16+16+8+4 ?= sizeof(struct CTABLE)

				// read all the actual shift values
				int a_size=ct.lim.lam*ct.lim.phi;

				ct.cvs=new LP[a_size];

				BinaryReader br=new BinaryReader(fid);
				for(int i=0; i<a_size; i++)
				{
					ct.cvs[i].lam=br.ReadSingle();
					ct.cvs[i].phi=br.ReadSingle();
				}
			}
			catch
			{
				ct.cvs=null;

#if DEBUG
				Console.Error.WriteLine("ctable loading failed on fread() - binary incompatible?");
#endif
				Proj.pj_ctx_set_errno(ctx, -38);
				return false;
			}

			return true;
		}

		//**********************************************************************
		//							nad_ctable_init()
		//
		//		Read the header portion of a "ctable" format grid.
		//**********************************************************************
		public static CTABLE nad_ctable_init(projCtx ctx, Stream fid)
		{
			try
			{
				CTABLE ct=new CTABLE();
				byte[] tmpID=new byte[80];

				// read the table header
				fid.Read(tmpID, 0, 80);

				BinaryReader br=new BinaryReader(fid);
				ct.ll.lam=br.ReadDouble();
				ct.ll.phi=br.ReadDouble();
				ct.del.lam=br.ReadDouble();
				ct.del.phi=br.ReadDouble();

				ct.lim.lam=br.ReadInt32();
				ct.lim.phi=br.ReadInt32();
				br.ReadInt32(); // FLP* cvs

				// do some minimal validation to ensure the structure isn't corrupt
				if(ct.lim.lam<1||ct.lim.lam>100000||ct.lim.phi<1||ct.lim.phi>100000)
				{
					Proj.pj_ctx_set_errno(ctx, -38);
					return null;
				}

				ct.id=Encoding.ASCII.GetString(tmpID);
				ct.id=ct.id.Trim(); // trim white space and newlines off id

				ct.cvs=null;

				return ct;
			}
			catch
			{
				Proj.pj_ctx_set_errno(ctx, -38);
				return null;
			}
		}

		//**********************************************************************
		//							nad_ctable2_load()
		//
		//		Load the data portion of a ctable2 formatted grid.
		//**********************************************************************
		public static bool nad_ctable2_load(projCtx ctx, CTABLE ct, Stream fid)
		{
			try
			{
				fid.Seek(160, SeekOrigin.Begin);

				// read all the actual shift values
				int a_size=ct.lim.lam*ct.lim.phi;

				ct.cvs=new LP[a_size];

				BinaryReader br=new BinaryReader(fid);

				if(IS_LSB)
				{
					for(int i=0; i<a_size; i++)
					{
						ct.cvs[i].lam=br.ReadSingle();
						ct.cvs[i].phi=br.ReadSingle();
					}
				}
				else
				{
					byte[] buf=br.ReadBytes(a_size*8);
					swap_words(buf, 0, 4, a_size*2);

					using(BinaryReader br2=new BinaryReader(new MemoryStream(buf)))
					{
						for(int i=0; i<a_size; i++)
						{
							ct.cvs[i].lam=br2.ReadSingle();
							ct.cvs[i].phi=br2.ReadSingle();
						}
					}
				}

				return true;
			}
			catch
			{
				ct.cvs=null;

#if DEBUG
				Proj.pj_log(ctx, PJ_LOG.ERROR, "ctable2 loading failed on fread() - binary incompatible?");
#endif

				Proj.pj_ctx_set_errno(ctx, -38);
				return false;
			}
		}

		//**********************************************************************
		//							nad_ctable2_init()
		//
		//		Read the header portion of a "ctable2" format grid.
		//**********************************************************************
		public static CTABLE nad_ctable2_init(projCtx ctx, Stream fid)
		{
			try
			{
				CTABLE ct=new CTABLE();
				byte[] header=new byte[160];

				// read the table header
				fid.Read(header, 0, 160);

				if(!IS_LSB)
				{
					swap_words(header, 96, 8, 4);
					swap_words(header, 128, 4, 2);
				}

				MemoryStream mem=new MemoryStream(header);
				BinaryReader br=new BinaryReader(mem);
				byte[] signature=new byte[16];
				byte[] tmpID=new byte[80];

				signature=br.ReadBytes(16);
				string sig=Encoding.ASCII.GetString(signature);
				if(sig.Substring(0, 9)!="CTABLE V2")
				{
					Proj.pj_log(ctx, PJ_LOG.ERROR, "ctable2 - wrong header!");
					Proj.pj_ctx_set_errno(ctx, -38);
					return null;
				}

				// read the table header
				tmpID=br.ReadBytes(80);

				ct.ll.lam=br.ReadDouble();
				ct.ll.phi=br.ReadDouble();
				ct.del.lam=br.ReadDouble();
				ct.del.phi=br.ReadDouble();

				ct.lim.lam=br.ReadInt32();
				ct.lim.phi=br.ReadInt32();

				// do some minimal validation to ensure the structure isn't corrupt
				if(ct.lim.lam<1||ct.lim.lam>100000||ct.lim.phi<1||ct.lim.phi>100000)
				{
					Proj.pj_ctx_set_errno(ctx, -38);
					return null;
				}

				ct.id=Encoding.ASCII.GetString(tmpID);
				ct.id=ct.id.Trim(); // trim white space and newlines off id

				ct.cvs=null;

				return ct;
			}
			catch
			{
				Proj.pj_ctx_set_errno(ctx, -38);
				return null;
			}
		}

		//**********************************************************************
		//								nad_init()
		//
		//		Read a datum shift file in any of the supported binary formats.
		//**********************************************************************
		public static CTABLE nad_init(projCtx ctx, string name)
		{
			ctx.last_errno=0;

			// --------------------------------------------------------------------
			//		Open the file using the usual search rules.
			// --------------------------------------------------------------------

			Stream fid=Proj.pj_open_lib(ctx, name, FileAccess.Read);
			if(fid==null) return null;

			CTABLE ct=nad_ctable_init(ctx, fid);
			if(ct!=null)
			{
				if(!nad_ctable_load(ctx, ct, fid))
				{
					nad_free(ct);
					ct=null;
				}
			}

			fid.Close();
			return ct;
		}

		//**********************************************************************
		//								nad_free()
		//
		//		Free a CTABLE grid shift structure produced by nad_init().
		//**********************************************************************
		public static void nad_free(CTABLE ct)
		{
			if(ct!=null) ct.cvs=null;
		}
	}
}
