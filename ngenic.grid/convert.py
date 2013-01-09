import ConfigParser
import argparse
import numpy
import os.path
import StringIO
parser = argparse.ArgumentParser(description="Assembling IC files from raw dumps")
parser.add_argument("paramfile", 
               help="the same paramfile fed to ngenic.grid.")
parser.add_argument("format", choices=['gadget', 'ramses'], 
               help="Format of the output. \
                     gadget : to write a gadget IC, \
                     ramses : to write a ramses IC, only the first zoom region is written due to ramses limitation")
parser.add_argument("-N", "--num-particles", 
               help="Number of DM particles per file. \
                     For levels with Gas the total number in a file will be twice of this",
                    default=1024 * 1024 * 20, dest='Npar', type=int);
parser.add_argument("-o", "--prefix", 
               help="Prefix of the IC files. \
                     Do not forget to put a filename base. \
                     gadget : IC.%%d, one file per zoom level. \
                     ramses : IC-%%d-%%s, several files per zoom level",
                                  default="IC")
parser.add_argument("-q", '--cubic', 
               help="Use an axis-aligned box for the zoom region, instead of spheres", 
                   action='store_true', default=False)
parser.add_argument('--no-displacement', dest='nodisp', 
                help="Do not displace the position of particles. \
                      Physics will be corrupted. Do not use", 
                      action='store_true', default=False)

def parseargs(parser):
  args = parser.parse_args()
  config = ConfigParser.ConfigParser()
  config.readfp(StringIO.StringIO(
"""[Cosmology]
h = 0.72
H = 0.1
G = 43007.1
C = 3e5
OmegaB = 0.044
OmegaM = 0.26
OmegaL = 0.74
[Units]
UnitLength_in_cm = 3.085678e21
UnitMass_in_g = 1.989e43
UnitVelocity_in_cm_per_s = 1e5
"""))
  str = file(args.paramfile).read().replace(';', ',').replace('#', ';')
  print str
  config.readfp(StringIO.StringIO(str))
  args.datadir = config.get("IO", "datadir")
  args.BoxSize = config.getfloat("IC", "BoxSize")
  args.a = config.getfloat("IC", "a")
  args.OmegaM = config.getfloat("Cosmology", "OmegaM")
  args.OmegaB = config.getfloat("Cosmology", "OmegaB")
  args.OmegaL = config.getfloat("Cosmology", "OmegaL")
  args.h = config.getfloat("Cosmology", "h")
  args.UnitLength_in_cm = config.getfloat("Unit", "UnitLength_in_cm")
  args.UnitVelocity_in_cm_per_s = config.getfloat("Unit", "UnitVelocity_in_cm_per_s")
  args.UnitMass_in_g = config.getfloat("Unit", "UnitMass_in_g")
  args.UnitTime_in_s = args.UnitLength_in_cm / args.UnitVelocity_in_cm_per_s;
  args.G = config.getfloat("Cosmology", "G")
  args.H = config.getfloat("Cosmology", "H")

  Regions = []
  for name, value in config.items("Regions"):
    all = [float(t) for t in value.replace(';', ',').split(',')]
    if len(all) == 4:
      x, y, z, s = all
      sx = s
      sy = s
      sz = s
    elif len(all) == 6:
      x, y, z, sx, sy, sz = all
  
    Regions.append(((x, y, z), (sx, sy, sz)))
  args.Regions = numpy.array(Regions, dtype=[('center', ('f8', 3)), ('size', ('f8', 3))])
  
  Levels = []
  
  args.Scale = {}
  args.Meta = {}
  args.DMptype = {}
  args.MakeGas = {}
  args.DownSample = {}
  for name, value in config.items("Levels"):
    sp = value.replace(';', ',').split(',')
    good = False
    if len(sp) < 2 or len(sp) > 4:
      raise 'paramfile format error, need 2/3/4 entries per Level [Levels]'

    if len(sp) >= 2:
      n, s = int(sp[0]), float(sp[1])
      # the default is 
      p = -1
      g = None
      d = 1
    if len(sp) >= 3:
      if 'g' in sp[2]:
        p = int(sp[2].replace('g', ''))
        g = True
      else:
        p = int(sp[2])
        g = False
    if len(sp) == 4:
      d = int(sp[3])

    Levels.append(n)
    args.Scale[n] = s
    args.DMptype[n] = p
    args.MakeGas[n] = g
    args.DownSample[n] = d
    args.Meta[n] = read_meta(n, args.datadir)
  args.Levels = numpy.array(Levels, dtype='i8')
  args.Levels.sort()
  for Nmesh in args.Levels:
    if args.DMptype[Nmesh] == -1:
      if Nmesh == args.Levels.max():
        args.DMptype[Nmesh] = 1
        args.MakeGas[Nmesh] = True
      else:
        args.DMptype[Nmesh] = 2
        args.MakeGas[Nmesh] = False
      
  return args

def read_meta(Nmesh, datadir):
    meta = {'Offset':{}, 'Size':{}}
    metalines = file('%s/meta-%d' % (datadir, Nmesh)).read()
    exec(metalines, meta)
    return meta

def writerecord(file, data):
    foo = numpy.empty(1, 'i4')
    foo[...] = data.nbytes
    foo.tofile(file)
    data.tofile(file)
    foo.tofile(file)


args = parseargs(parser)

def setup_level(Nmesh):
  if args.Scale[Nmesh] == 0: 
    return setup_baselevel(Nmesh)
  else:
    return setup_innerlevel(Nmesh)
  
def readblock(Nmesh, block, dtype):
  header = {}
  exec(
    file('%s/%s-%d.header' % (args.datadir, block, Nmesh)).read(), 
    header)
  NTask = header['NTask']
  DownSample = header['DownSample']
  assert DownSample == args.DownSample[Nmesh]

  width1 = len(str(DownSample))
  width2 = len(str(NTask))
  content = []

  f = 0
  for i in range(DownSample):
    for fid in range(NTask):
      if DownSample > 1:
        fname = '%s/%s-%d.%0*d-%0*d' %(
                 args.datadir, 
                 block, Nmesh, width1, i, width2, fid)
      else:
        fname = '%s/%s-%d.%0*d' %(
                 args.datadir, 
                   block, Nmesh, width2, fid)
      if os.path.exists(fname):
        content.append(numpy.fromfile(fname, dtype))
        f = f + 1
      #skip the tasks that do not dump any files
  print Nmesh, block, f, 'files read'
  return numpy.concatenate(content)

def setup_innerlevel(Nmesh):
  index = readblock(Nmesh, 'index', 'i8')
  r = readblock(Nmesh, 'region', 'u4')
  dispx = readblock(Nmesh, 'dispx', 'f4')
  dispy = readblock(Nmesh, 'dispy', 'f4')
  dispz = readblock(Nmesh, 'dispz', 'f4')
  delta = readblock(Nmesh, 'delta', 'f4')
  meta = args.Meta[Nmesh]
  assert len(index) == len(r)
  assert len(index) == len(dispx)
  assert len(index) == len(dispy)
  assert len(index) == len(dispz)
  assert len(index) == len(delta)
  alldata = []
  for i, region in enumerate(args.Regions):
    mask = r == i
    data = numpy.empty(mask.sum(), dtype=[
           ('gps', ('i8', 3)), 
           ('disp', ('f4', 3)),
           ('delta', 'f4')])
    data['gps'][:, 0], data['gps'][:, 1], data['gps'][:, 2] = numpy.unravel_index(index[mask], meta['Size'][i])
    data['gps'] += meta['Offset'][i]
    data['disp'][:, 0] = dispx[mask]
    data['disp'][:, 1] = dispy[mask]
    data['disp'][:, 2] = dispz[mask]
    data['delta'] = delta[mask]
    alldata.append(data)
  data = numpy.concatenate(alldata)
  # now gonna remove the overlaps and make sure points are within the mesh.
  numpy.remainder(data['gps'], Nmesh, data['gps'])
  index = numpy.ravel_multi_index(data['gps'].T, (Nmesh, Nmesh, Nmesh))
  u, ui = numpy.unique(index, return_index=True)
  return data[ui]

def setup_baselevel(Nmesh):
  data = numpy.empty((Nmesh, Nmesh, Nmesh), dtype=[
           ('gps', ('i8', 3)), 
           ('disp', ('f4', 3)),
           ('delta', 'f4')])

  grid1d = numpy.arange(Nmesh);

  data['gps'][..., 2] = grid1d[None, None, :]
  data['gps'][..., 1] = grid1d[None, :, None]
  data['gps'][..., 0] = grid1d[:, None, None]

  data = data.reshape(-1)
  for i, block in enumerate(['dispx', 'dispy', 'dispz']):
    content = readblock(Nmesh, block, 'f4')
    data['disp'][..., i] = content

  data['delta'] = readblock(Nmesh, 'delta', 'f4')

  return data

def ramses():
  def write_ramses(header, Nmesh, name, buffer):
    with file('%s%d-%s' % (args.prefix, Nmesh, name), 'w') as icfile:
      buffer = buffer.reshape(header['np'])
      writerecord(icfile, header)
      for row in buffer.transpose((2, 0, 1)):
        writerecord(icfile, row.ravel('F'))

  RHEADER = numpy.dtype([
   ('np', ('i4', 3)), ('dx', 'f4'), 
   ('xo', ('f4', 3)), ('astart', 'f4'), 
   ('omegam', 'f4'), ('omegav', 'f4'), ('h0', 'f4')])
  fac = 1.0 / args.h * args.UnitLength_in_cm / 3.08567758e24
  boxsize = args.BoxSize * fac
  for i, Nmesh in enumerate(args.Levels):
    header = numpy.zeros(None, RHEADER)
    if args.Scale[Nmesh] == 0:
      header['np'][:] = Nmesh
      header['xo'][:] = 0
      mask = None
    else:
      header['np'][:] = args.Meta[Nmesh]['Size'][0]
      header['xo'][:] = header['dx'] * args.Meta[Nmesh]['Offset'][0]
      r = readblock(Nmesh, 'region', 'u1')
      mask = r == 0
    header['dx'] = boxsize / Nmesh
    header['astart'] = args.a
    header['omegam'] = args.OmegaM
    header['omegav'] = args.OmegaL
    header['h0'] = args.h * 100.0
    print Nmesh, header

    delta = readblock(Nmesh, 'delta', 'f4')
    if mask is not None: delta = delta[mask]
    write_ramses(header, Nmesh, 'deltab', delta)
    delta = None
    for i, block in enumerate(['dispx', 'dispy', 'dispz']):
      disp = readblock(Nmesh, block, 'f4')
      if mask is not None: disp = disp[mask]
      disp *= fac
      write_ramses(header, Nmesh, block, disp)
      disp /= args.Meta[Nmesh]['vfact']
      write_ramses(header, Nmesh, ['velcx', 'velcy', 'velcz'][i], disp)

def gadget():
  GHEADER = numpy.dtype([
      ('N', ('u4', 6)),
      ('mass', ('f8', 6)),
      ('time', 'f8'),
      ('redshift', 'f8'),
      ('flag_sfr', 'i4'),     # unused 
      ('flag_feedback', 'i4'),  # unused
      ('Ntot_low', ('u4', 6)),
      ('flag_cool', 'i4'),  # unused
      ('Nfiles', 'i4'),
      ('boxsize', 'f8'),
      ('OmegaM', 'f8'),
      ('OmegaL', 'f8'),
      ('h', 'f8'),
      ('flag_sft', 'i4'),
      ('flag_met', 'i4'),
      ('Ntot_high', ('u4', 6)),
      ('flag_entropy', 'i4'),
      ('flag_double', 'i4'),
      ('flag_ic_info', 'i4'),
      ('flag_lpt_scalingfactor', 'i4'),
      ('unused', ('u4', 12)),
  ])

  def dig(data, Nmesh, scale, Align=None):
    """ dig a hole on Nmesh points data, for NmeshDig
        returns the mask of the hole,
        align to Align, if omitted use Nmesh
    """
    if Align is None: Align = Nmesh
    sep = args.BoxSize / Align
    assert Nmesh % Align == 0
    mask2 = numpy.zeros(len(data), dtype='?')
    for region in args.Regions:
      gridcenter = region['center'] / sep
      dis = data['gps'] // (Nmesh // Align) - gridcenter
      # use the if any part of the mesh overlaps the sphere,
      dis[dis < - 0.5] += 1
      dis *= sep
      dis += args.BoxSize * 0.5
      numpy.remainder(dis, args.BoxSize, dis)
      dis -= args.BoxSize * 0.5
      dis /= (0.5 * region['size'] * scale)
      if not args.cubic:
        dis **= 2
        dis = dis.sum(axis=-1)
        mask2 |= (dis < 1.0)
      else:
        numpy.abs(dis, dis)
        mask2 |= (dis < 1.0).all(axis=-1)
    return mask2

  def rewrite_header(fname, header):
    with file(fname, 'r+') as icfile:
      writerecord(icfile, header)
    
  def write_gadget(fname, data):
    Nmesh = args.Levels[i]
    ptype = args.DMptype[Nmesh]
    makegas = args.MakeGas[Nmesh]

    with file(fname, 'w') as icfile:
      critical_mass = 3 * args.H ** 2 / (8 * numpy.pi * args.G) * (args.BoxSize / Nmesh) ** 3
  
      header = numpy.zeros(None, GHEADER)
      header['time'] = args.a
      header['redshift'] = 1. / args.a - 1
      header['boxsize'] = args.BoxSize
      header['OmegaM'] = args.OmegaM
      header['OmegaL'] = args.OmegaL
      header['h'] = args.h
      Npar = len(data)
  
      header['N'][ptype] = Npar
      usemasstab = numpy.sum([args.DMptype[n] == ptype for n in args.Levels]) == 1
      if args.MakeGas[Nmesh]:
        header['N'][0] = Npar
        if usemasstab:
          header['mass'][ptype] = (args.OmegaM - args.OmegaB) * critical_mass
      else:
        if usemasstab:
          header['mass'][ptype] = args.OmegaM * critical_mass
  
      writerecord(icfile, header)
  
      # position
      pos = numpy.float32(data['gps']) * args.BoxSize
      pos /= Nmesh
      if not args.nodisp:
        pos += data['disp']
      if makegas: 
        pos = numpy.tile(pos, (2, 1))
        # gas
        pos[:Npar] -= 0.5 * (1 - args.OmegaB / args.OmegaM) * args.BoxSize / Nmesh
        # dm
        pos[Npar:] += 0.5 * args.OmegaB / args.OmegaM * args.BoxSize / Nmesh
  
      numpy.remainder(pos, args.BoxSize, pos)
  
      writerecord(icfile, pos)
      pos = None
  
      # velocity
      vel = data['disp'] * args.Meta[Nmesh]['vfact']
      if makegas: 
        vel = numpy.tile(vel, (2, 1))
      writerecord(icfile, vel)
      vel = None
  
      # id
      if makegas:
        id = numpy.empty(Npar + Npar, dtype=('i4', 2))
        id[:, 0] = numpy.arange(Npar + Npar)
      else:
        id = numpy.empty(Npar, dtype=('i4', 2))
        id[:, 0] = numpy.arange(Npar)
      id[:, 1] = len(args.Levels) - i - 1
      writerecord(icfile, id)
      id = None
  
      # mass
      if usemasstab:
        if makegas:
          mass = numpy.empty(Npar, dtype='f4')
          mass[:] = args.OmegaB * critical_mass
          writerecord(icfile, mass)
          mass = None 
      else:
        if makegas:
          mass = numpy.empty(Npar + Npar, dtype='f4')
          mass[:Npar] = args.OmegaB * critical_mass
          mass[Npar:] = (args.OmegaM - args.OmegaB) * critical_mass
        else:
          mass = numpy.empty(Npar, dtype='f4')
          mass[:Npar] = args.OmegaM * critical_mass
        writerecord(icfile, mass)
        mass = None 
      # ie ( all zeros)
      if makegas:
        ie = numpy.zeros(Npar, dtype='f4')
        writerecord(icfile, ie)
        ie = None
    return header

  H = []
  fid = 0
  for i, Nmesh in enumerate(args.Levels):
    print i, Nmesh
    data = setup_level(Nmesh)
    if i < len(args.Levels) - 1:
      # if there is a refine level, dig a hole from the current level
      # leaving space for the next level.
      NmeshNext = args.Levels[i + 1]
      assert NmeshNext % Nmesh == 0
      holemask = dig(data, Nmesh, scale=args.Scale[NmeshNext])
      if args.DownSample[NmeshNext] > 1:
        # downsample levels need to carry on the interpolation
        # of previous level
        # we use the corresponding upper level values.
        # It really should have been a fourier interploation.
        # first save them to a temp place
        # then load them after the 'fill' step is done.
        dataPrevTemp = data[holemask] 
        indexPrevTemp = numpy.ravel_multi_index(dataPrevTemp['gps'].T, (Nmesh, Nmesh, Nmesh))
        arg = indexPrevTemp.argsort()
        dataPrevTemp = dataPrevTemp[arg]
        indexPrevTemp = indexPrevTemp[arg]
        print Nmesh, 'hole', dataPrevTemp['gps'].min(axis=0), dataPrevTemp['gps'].max(axis=0), indexPrevTemp.min(), indexPrevTemp.max(), indexPrevTemp
        
        arg = None
      data = data[~holemask]
    if i > 0:
      # if there is a parent level, clip the current level
      # to fill in the space of the parent level
      NmeshPrev = args.Levels[i - 1]
      assert Nmesh % NmeshPrev == 0
      fillmask = dig(data, Nmesh, scale=args.Scale[Nmesh], Align=NmeshPrev)
      print Nmesh, 'before fill', data['gps'].min(axis=0), data['gps'].max(axis=0)
      data = data[fillmask]
      print Nmesh, 'after fill', data['gps'].min(axis=0), data['gps'].max(axis=0)
      if args.DownSample[Nmesh] > 1:
        # the DownSample interpolation.
        # first find the corresponding lowerlevel points
        index = numpy.ravel_multi_index(data['gps'].T / (Nmesh / NmeshPrev),
                   (NmeshPrev, NmeshPrev, NmeshPrev))
        print Nmesh, 'index', index.min(), index.max()
        arg = indexPrev.searchsorted(index, 'left')
        if not (indexPrev[arg] == index).all():
#          print data[indexPrev[arg] != index]
          print indexPrev, index
          print 'assertion failure for these points', Nmesh, NmeshPrev
          raise Exception("")
        # then add them to the current level
        data['disp'] += dataPrev['disp'][arg]
        data['delta'] += dataPrev['delta'][arg]
        arg = None
        index = None
    if i < len(args.Levels) - 1:
      if args.DownSample[NmeshNext] > 1:
        indexPrev = indexPrevTemp
        dataPrev = dataPrevTemp

    nchunks = int(numpy.ceil(len(data) / (1.0 * args.Npar)))
    datasp = numpy.array_split(data, nchunks)
    for k in range(nchunks):
      fname = '%s.%d' % (args.prefix, fid + k)
      H.append((fname, write_gadget(fname, datasp[k])))
    fid = fid + nchunks

  Ntot = numpy.zeros(6, dtype='i8')
  masstab = numpy.zeros(6, dtype='f8')

  for fname, header in H:
    Ntot += header['N']
    mask = header['mass'] != 0
    masstab[mask] = header['mass'][mask]
  print 'total particles written', Ntot
  for fname, header in H:
    header['mass'] = masstab
    header['Ntot_low'][:] = (Ntot & 0xffffffff)
    header['Ntot_high'][:] = (Ntot >> 32)
    header['Nfiles'] = len(H)
  
    print fname, header['N']
    rewrite_header(fname, header)

def main():
  if args.format == 'gadget':
    gadget()
  elif args.format == 'ramses':
    ramses()

main()
