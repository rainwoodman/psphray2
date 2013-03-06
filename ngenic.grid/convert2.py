import numpy
from scipy.ndimage import map_coordinates
# this requires numpy > 1.6.0
from args import parseargs
from myio import readchunk, writechunk, yieldfilename, write_gadget_one, rewrite_header

def interpolate(h, hnext):
  """ interpolate h into hnext """
  # round down
  exbot = (hnext['Offset'] * h['Nmesh']) // hnext['Nmesh']
  # round up
  extop = ((hnext['Offset'] + hnext['Size']) * h['Nmesh'] + h['Nmesh'] - 1)// hnext['Nmesh']

  print h['Nmesh'], 'box', exbot, extop

  overlay_list = []
  for fn in list(yieldfilename(h)):
    chunk = readchunk(h, fn, exbot, extop)
    ipos = chunk['ipos']
    print 'ipos stat', 'min', ipos.min(axis=0), 'max', ipos.max(axis=0)
    overlay_list.append(chunk)
  overlay = numpy.concatenate(overlay_list)
  overlay_list = []
  index = numpy.ravel_multi_index((overlay['ipos'] - exbot).T, extop - exbot)
  overlay = overlay[index.argsort()]

  for fn in list(yieldfilename(hnext)):
    # read in the current level displacement
    chunk = readchunk(hnext, fn, None, None, extra='1')
    ipos = chunk['ipos']
    ipos = 1.0 * ipos * h['Nmesh'] / hnext['Nmesh'] - exbot
    src = overlay['delta'].reshape(*(extop-exbot))
    chunk['delta'] += map_coordinates(src, ipos.T)
    for d in range(3):
      src = overlay['disp'][:, d].reshape(*(extop-exbot))
      chunk['disp'][:, d] += map_coordinates(src, ipos.T)
    # write to the full displacement
    writechunk(hnext, fn, chunk)

def main():
  A = parseargs()
  if A.action == 'interpolate':
    main_interpolate(A)
  if A.action == 'gadget':
    main_gadget(A)

def main_interpolate(A):
  Nmesh = A.L['Nmesh']
  for i in range(1, len(A.L)):
    h = A.META[Nmesh[i - 1]]
    hnext = A.META[Nmesh[i]]
    if hnext['DownSample'] > 1:
      interpolate(h, hnext)

def main_gadget(A):
  Nmesh = A.L['Nmesh']
  fid = 0
  F = []
  for i in range(len(A.L)):
    h = A.META[Nmesh[i]]

    for fn in list(yieldfilename(h)):
      chunk = readchunk(h, fn, None, None)
      ipos = chunk['ipos']
      if i != 0:
        IncludeCenter = h['ZoomCenter']
        IncludeRadius = h['ZoomRadius']
        includemask = selectmask(chunk, 
               A.BoxSize, 
               Nmesh[i],
               Nmesh[i - 1],
               IncludeCenter, IncludeRadius, A.cubic)
        chunk = chunk[includemask]
      if i != len(A.L) - 1:
        hnext = A.META[Nmesh[i + 1]]
        ExcludeCenter = hnext['ZoomCenter']
        ExcludeRadius = hnext['ZoomRadius']

        excludemask = selectmask(chunk,
               A.BoxSize,
               Nmesh[i], 
               Nmesh[i],
               ExcludeCenter, ExcludeRadius, A.cubic)
        chunk = chunk[~excludemask]
      print 'ipos stat', 'min', ipos.min(axis=0), 'max', ipos.max(axis=0)
      F.append(write_gadget_one('%s.%d' % (A.prefix, fid), h, chunk))
      fid = fid + 1

  Ntot = numpy.zeros(6, dtype='i8')

  Ntot[:] = numpy.sum([ header['N'] for fname, header in F], axis=0)
  print 'total particles written', Ntot

  for fname, header in F:
    header['mass'] = 0
    header['Ntot_low'][:] = (Ntot & 0xffffffff)
    header['Ntot_high'][:] = (Ntot >> 32)
    header['Nfiles'] = len(F)
  
    print 'touching up',  fname
    rewrite_header(fname, header)

def selectmask(chunk, BoxSize, Nmesh, Nalign, center, radius, cubic):
  assert Nmesh % Nalign == 0
  AlignSpacing = BoxSize / Nalign
  gridcenter = center / AlignSpacing
  halfsize = 1. / radius
  dis = chunk['ipos'] // (Nmesh // Nalign) - gridcenter
  dis[dis < -0.5] += 1
  dis *= AlignSpacing
  dis += BoxSize * 0.5
  numpy.remainder(dis, BoxSize, dis)
  dis -= BoxSize * 0.5
  dis *= halfsize
  if not cubic:
    dis **= 2
    dis = dis.sum(axis=-1)
    mask = dis < 1.0
  else:
    numpy.abs(dis, dis)
    mask = (dis < 1.0).all(axis=-1)
  return mask
main()
