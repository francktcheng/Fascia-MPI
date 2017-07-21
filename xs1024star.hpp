

struct xs1024star_t {
  unsigned long s[16];
  long p;
} ;

unsigned long xs1024star_next(xs1024star_t* xs) 
{
   const unsigned long s0 = xs->s[xs->p];
   unsigned long s1 = xs->s[xs->p = (xs->p + 1) & 15];
   s1 ^= s1 << 31;
   xs->s[xs->p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30);
   return xs->s[xs->p] * (unsigned long)(1181783497276652981U);
}

double xs1024star_next_real(xs1024star_t* xs) 
{
   const unsigned long s0 = xs->s[xs->p];
   unsigned long s1 = xs->s[xs->p = (xs->p + 1) & 15];
   s1 ^= s1 << 31;
   xs->s[xs->p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30);
   double ret = (double)(xs->s[xs->p] * (unsigned long)(1181783497276652981U));   
   return ret /= (double)((unsigned long)(18446744073709551615U));
}

void xs1024star_seed(unsigned long seed, xs1024star_t* xs) 
{
  for (unsigned long i = 0; i < 16; ++i)
  {
    unsigned long z = (seed += (unsigned long)(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * (unsigned long)(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * (unsigned long)(0x94D049BB133111EB);
    xs->s[i] = z ^ (z >> 31);
  }
  xs->p = 0;
}
