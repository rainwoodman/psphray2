import numpy
from zoom import zoom

def trilinearinterp(x, xs, ys, dims, mode='wrap', out=None):
  """ interpolate x from grid points (xs, ys),
      the grid points are on a mesh of dims.
      trilinear.
      x can be non-integer.
      xs must be integer.
      out type will be in ys type, or f8.
      internally the calculation is always done in 'f8'
      
      x is subdivied to chunks of 8192 elements, so that this won't use too much memory.
      will raise ValueError if the input mesh do not cover x.
      mode can be 'raise', 'clip' or 'wrap', deciding the boundary condition.
  """
  
  if out is None: 
    if hasattr(ys, 'dtype'): dtype = ys.dtype
    else: dtype = 'f8'
    if hasattr(ys, 'shape') and len(ys.shape) > 1:
      shape = [len(x)] + list(ys.shape[1:])
    else:
      shape = len(x)
    
    out = numpy.empty(shape, dtype)
  if len(x) == 0: return out
  vec = numpy.array([
               [0, 0, 0],
               [0, 0, 1],
               [0, 1, 0],
               [1, 0, 0],
               [0, 1, 1],
               [1, 0, 1],
               [1, 1, 0],
               [1, 1, 1]], dtype='i8')
  sorted_raveled_xs = numpy.ravel_multi_index(numpy.asarray(xs).T, dims, mode=mode)
  arg = sorted_raveled_xs.argsort()
  ys = ys[arg]
  sorted_raveled_xs = sorted_raveled_xs[arg]
  del arg

  def loop(x, out):
    corner = numpy.empty(len(x), dtype=('i4', 3))
    numpy.floor(x, corner)
    residual = x - corner

    tmp = numpy.zeros(out.shape, dtype='f8')
    for v in vec:
      weight = ((2 * v - 1) * residual + (1 - v)).prod(axis=-1)
      corner += v
      raveled_x = numpy.ravel_multi_index(corner.T, dims, mode=mode)
      corner -= v
      ind = sorted_raveled_xs.searchsorted(raveled_x)
      found = sorted_raveled_xs.take(ind, mode='clip') == raveled_x
      if not found.all():
        raise ValueError("Some points from x is not covered in xs")
      tmp += (ys[ind].T * weight).T
    
    out[:] = tmp
  for i in range(0, len(x), 65536):
    loop(x[i:i+65536], out[i:i+65536])

  return out
def upscaleby2(src):
  return zoom(src, 2, order=5, mode='wrap')
  newshape = [ s * 2 for s in src.shape]
  x = numpy.array(numpy.unravel_index(range(len(src.flat) * 8),
       newshape)).T
  xs = numpy.array(numpy.unravel_index(range(len(src.flat)),
       src.shape)).T

  return trilinearinterp(x * 0.5, xs, src.ravel(), dims=src.shape, mode='wrap').reshape(newshape)

def power(delta, K0, Dplus):
  delta_k = numpy.fft.rfftn(delta) / numpy.prod(numpy.float64(delta.shape))
  i, j, k = numpy.unravel_index(range(len(delta_k.flat)), delta_k.shape)
  i = delta_k.shape[0] // 2 - numpy.abs(delta_k.shape[0] // 2 - i)
  j = delta_k.shape[1] // 2 - numpy.abs(delta_k.shape[1] // 2 - j)
  k = (i ** 2 + j ** 2 + k ** 2) ** 0.5
  dig = numpy.digitize(k, numpy.arange(int(k.max())))
  kpk = k.reshape(delta_k.shape) * (numpy.abs(delta_k) ** 2) * K0 ** -3 * Dplus ** 2
  return numpy.arange(dig.max()+1) * K0, numpy.bincount(dig, weights=kpk.ravel()) / numpy.bincount(dig, weights=k.ravel())
  
def main():  
  dims = (10, 10, 10)
  xs = numpy.array(numpy.unravel_index(numpy.arange(10 ** 3), dims)).T
  ys = 1.0 * numpy.mod(xs[:, 0], 10)

  x = [[i, 0, i] for i  in numpy.linspace(0, 20, 41)]
  y = trilinearinterp(x, xs, ys, dims)

  print zip(x, y)

if __name__ == '__main__':
  # test
  main()

