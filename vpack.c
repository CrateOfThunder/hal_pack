/* ------------------------------------------------------------------ */
/*                                                                    */
/*                   HISTORICAL PRESERVATION OBJECT                   */
/*                                                                    */
/* ------------------------------------------------------------------ */





/*
 * This software is available to the public as a historical object.
 * It must not be used beyond personal/educational application(s),
 * and, subsequently, any premeditated intents to sell or incorporate
 * the sofware within other projects contradicting the aforementioned
 * conditions will be considered malignant & condemnable.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE & NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/* ------------------------------------------------------------------ */
/* Vpack - CGB Pocket Monsters Revision                               */
/*                    vpack.c -- gfx [de]compressor                   */
/* Rev.                        28JAN2025               CrateOfThunder */
/* Ver. 1.XX                   ??SEP1998    S. Iwata / HAL Laboratory */
/* Ver. 1.60                   13DEC1994    S. Iwata / HAL Laboratory */
/* Ver. 1.00                   ??JUN1989    S. Iwata / HAL Laboratory */
/* ------------------------------------------------------------------ */
/*                                                                    */
/*       ANY PERCEIVED PERSONAL CLAIM OF INTELLECTUAL PROPERTY        */
/*         NEITHER IMPLIED NOR INFERRABLE PER THE ACCOUNT OF          */
/*                           CrateOfThunder                           */
/*                                                                    */
/* ------------------------------------------------------------------ */

/*
 * gcc -x c -std=c99 -Wall -Wextra -Wpedantic -Werror -Os -s
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <stdint.h>

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;

/* ------------------------------------------------------------------ */
/*                                                                    */
/*                  GENERIC C VECTOR IMPLEMENTATION                   */
/*                                                                    */
/* ------------------------------------------------------------------ */

/* Define the dynamic array structure */
typedef struct {
  void *d;   /* data */
  size_t s;  /* size */
  size_t c;  /* capacity */
  size_t es; /* elem_size */
} Vec;

/* Initialize vector */
void v_init(Vec *v, size_t es)
{
  v->s = 0;
  v->c = 4;
  v->es = es;
  v->d = malloc(es * v->c);

  if (v->d == NULL) {
    fprintf(stderr, "v_init: alloc error\n");
    exit(EXIT_FAILURE);
  }

  return;
}

/* Free vector memory */
void v_free(Vec *v)
{
  free(v->d);
  v->d = NULL;
  v->s = 0;
  v->c = 0;
  v->es = 0;
  return;
}

/* Resize vector */
static void v_resize(Vec *v, size_t nc)
{
  void *nd = realloc(v->d, v->es * nc);

  if (nd == NULL) {
    fprintf(stderr, "v_resize: alloc error\n");
    exit(EXIT_FAILURE);
  }

  v->d = nd;
  v->c = nc;
  return;
}

/* Append element */
void v_push(Vec *v, const void *val)
{
  if (v->s >= v->c) { v_resize(v, v->c * 2); }

  memcpy((char *)v->d + v->s * v->es, val, v->es);
  v->s++;
  return;
}

/* Get element */
void *v_get(Vec *v, size_t i)
{
  if (i >= v->s) {
    fprintf(stderr, "v_get: index error\n");
    exit(EXIT_FAILURE);
  }

  return (char *)v->d + i * v->es;
}

/* Set element */
void v_set(Vec *v, size_t i, const void *val)
{
  if (i >= v->s) {
    fprintf(stderr, "v_set: index error\n");
    exit(EXIT_FAILURE);
  }

  memcpy((char *)v->d + i * v->es, val, v->es);
  return;
}

/* Insert element */
void v_ins(Vec *v, size_t i, const void *val)
{
  if (i > v->s) {
    fprintf(stderr, "v_ins: index error\n");
    exit(EXIT_FAILURE);
  }

  if (v->s >= v->c) { v_resize(v, v->c * 2); }

  memmove((char *)v->d + (i + 1) * v->es, 
          (char *)v->d + i * v->es, 
          (v->s - i) * v->es);
  memcpy((char *)v->d + i * v->es, val, v->es);
  v->s++;
  return;
}

/* Erase element */
void v_erase(Vec *v, size_t i)
{
  if (i >= v->s) {
    fprintf(stderr, "v_erase: index error\n");
    exit(EXIT_FAILURE);
  }

  memmove((char *)v->d + i * v->es, 
          (char *)v->d + (i + 1) * v->es, 
          (v->s - i - 1) * v->es);
  v->s--;
  return;
}

/* Get vector size */
size_t v_size(Vec *v) { return v->s; }

/* Get vector capacity */
size_t v_cap(Vec *v) { return v->c; }

/* Check if vector is empty */
int v_empty(Vec *v) { return v->s == 0; }

/* Clear vector */
void v_clear(Vec *v) { v->s = 0;  return; }



/* ------------------------------------------------------------------ */
/*                                                                    */
/*                       REFACTORED SOURCE CODE                       */
/*                                                                    */
/* ------------------------------------------------------------------ */
#define CGB_VRAM_SIZE 128*3*2*16
#define DATASIZELIMIT 65536      /* Max DAT that can be handled: 64KB */
#define EMPTY (u16)DATASIZELIMIT /* Marks a space in the hash table */
#define HASHTABLESIZE 131111     /* Prime numbers: 524287 512KB */
#define LONG_ENOUGH 64
#define NORMAL     (u8)0   /* Uncompressed normal data */
#define RUNLENGTH  (u8)1   /* 1 data sequence */
#define RUNLENGTH2 (u8)2   /* 2 data alternating continuous */
#define INCRUNLEN  (u8)3   /* Incremental Data Continuous */
#define REFERENCE  (u8)4   /* Normal Data Reference */
#define REFERLR    (u8)5   /* Left-right inversion data reference */
#define REFERUD    (u8)6   /* Refer to the upside-down data */
#define SPECIAL    (u8)7   /* Special Data */
#define ENDMARK    (u8)255 /* Data Termination */

static u32 m_algoPacked[7] = {1, 2, 3, 2, 3, 3, 3};
static u32 m_algoPackedZ[7] = {1, 2, 3, 1, 2, 2, 2};

/* PUBLIC */
Vec m_aBinData;
u32 m_CgxDataSize;
u8 m_aCgxData[CGB_VRAM_SIZE];
/* PRIVATE */
static u8 m_invertTable[256];           /* Left-right inversion table */
static u32 m_nrefHash[HASHTABLESIZE];
static u32 m_vrefHash[HASHTABLESIZE];
static u8 m_algorithm[DATASIZELIMIT];   /* Compression Algorithms */
static u32 m_compLength[DATASIZELIMIT]; /* Data LEN to be compressed */
static u32 m_reference[DATASIZELIMIT];  /* Reference Index */
static u32 m_algoPacked[7];
static u32 m_algoPackedZ[7];
static u32 m_packedSize;                /* Data LEN after compression */

static void outData(void)
{
  u32 i, j, k, len, loc;
  u8 tmp;

  v_clear(&m_aBinData);
  i = 0;

  while (i < m_CgxDataSize) {
    len = m_compLength[i];
    loc = i;

    do {
      j = (len > 1024) ? 1024 : len;

      if (j > 32) {
        /* Compression length 33 or more */
        tmp = (224 + (m_algorithm[i] << 2) + ((j - 1) >> 8));
        v_push(&m_aBinData, &tmp);
        tmp = ((j - 1) & 255);
        v_push(&m_aBinData, &tmp);
      }
      else {
        /* Compression length 1 to 32 */
        tmp = ((m_algorithm[i] << 5) + j - 1);
        v_push(&m_aBinData, &tmp);
      }

      switch (m_algorithm[i]) {
        case NORMAL:
          for (k = 0; k < j; k++) {
            v_push(&m_aBinData, &m_aCgxData[loc + k]);
          }

          loc += j;
          break;
        case RUNLENGTH:
          v_push(&m_aBinData, &m_aCgxData[i]);
          break;
        case INCRUNLEN:
          break;
        case RUNLENGTH2:
          v_push(&m_aBinData, &m_aCgxData[i]);
          v_push(&m_aBinData, &m_aCgxData[i + 1]);
          break;
        case REFERENCE:
        case REFERLR:
          /*if (i - m_reference[i] <= 128) {
            tmp = ((i - m_reference[i] - 1) | 128);
            v_push(&m_aBinData, &tmp);
          }
          else {
            tmp = (m_reference[i] >> 8);
            v_push(&m_aBinData, &tmp);
            tmp = (m_reference[i] & 255);
            v_push(&m_aBinData, &tmp);
          }

          break;*/
        case REFERUD:
          if (i - m_reference[i] <= 128) {
            tmp = ((i - m_reference[i] - 1) | 128);
            v_push(&m_aBinData, &tmp);
          }
          else {
            tmp = (m_reference[i] >> 8);
            v_push(&m_aBinData, &tmp);
            tmp = (m_reference[i] & 255);
            v_push(&m_aBinData, &tmp);
          }

          break;
      }

      len -= j;
    } while (len);

    i += m_compLength[i];
  }

  tmp = ENDMARK;
  v_push(&m_aBinData, &tmp); /* Write the end mark */
  return;
}

static void judge(void)
{
  u32 i, j;

  /* If compressing the DAT will be detrimental, it will be excluded. */
  for (i = 0; i < m_CgxDataSize; i++) {
    if (m_algorithm[i]) {
      for (j = 1; j < m_compLength[i]; j++) {
        if (m_compLength[i + j] > m_compLength[i]) {
          m_compLength[i] = j;

          if ((m_algorithm[i] == RUNLENGTH2) && (j & 1)) {
            m_compLength[i]--;
          }
        }
      }
    }
  }

  /* If data does not become shorter by simply compressing it,
     it is excluded from compression. */
  for (i = 0; i < m_CgxDataSize; i++) {
    if (m_compLength[i] <= m_algoPackedZ[m_algorithm[i]]) {
      m_algorithm[i] = NORMAL;
      m_compLength[i] = 0;
    }

    if ((m_algorithm[i] >= REFERENCE) && ((i - m_reference[i]) > 128)) {
      if (m_compLength[i] <= m_algoPacked[m_algorithm[i]]) {
        m_algorithm[i] = NORMAL;
        m_compLength[i] = 0;
      }
    }
  }

  /* Check for continuity of uncompressed data */
  for (i = 0; i < m_CgxDataSize; /**/) {
    if (m_algorithm[i] == NORMAL) {
      for (j = 1; (i + j) < m_CgxDataSize; j++) {
        if (m_algorithm[i + j] != NORMAL) {
          break;
        }
      }

      m_compLength[i] = j;
      i += j;
    }
    else {
      i += m_compLength[i];
    }
  }

  return;
}

/* Hash Functions */
static u32 hashf(u8 d0, u8 d1, u8 d2, u8 d3)
{
  u32 x = (u32)d0 + ((u32)d1 << 8) + ((u32)d2 << 16) + ((u32)d3 << 24);

  return (x % (u32)HASHTABLESIZE);
}

/* Register in a hash table */
static void registVHash(u32 x)
{
  u32 y, d = 1;

  y = (x >= 4) ? hashf(m_aCgxData[x - 0], m_aCgxData[x - 1],
                       m_aCgxData[x - 2], m_aCgxData[x - 3]) : 0;

  while (m_vrefHash[y] != EMPTY) {
    d += 2;

    if ((y += d) >= HASHTABLESIZE) {
      y -= HASHTABLESIZE;
    }
  }

  m_vrefHash[y] = x;
  return;
}

static void registNHash(u32 x, u32 y)
{
  u32 d = 1;

  while (m_nrefHash[y] != EMPTY) {
    d += 2;

    if ((y += d) >= HASHTABLESIZE) {
      y -= HASHTABLESIZE;
    }
  }

  m_nrefHash[y] = x;
  return;
}

/* Return values ​​sequentially from a hash table with open addressing */
#define NREF 0
#define VREF 1
static u32 scanHash(u32 h, int ref)
{
  static u32 lastH = HASHTABLESIZE, y, d;
  u32 r;

  if (h != lastH) {
    y = h;
    d = 1;
  }

  r = (ref == NREF) ? m_nrefHash[y] : m_vrefHash[y];
  lastH = (r == EMPTY) ? HASHTABLESIZE : h;
  d += 2;

  if ((y += d) >= HASHTABLESIZE) {
    y -= HASHTABLESIZE;
  }

  return r;
}

static void checkData(void)
{
  u32 i, j, k, l, r, hn, hh;
  u8 *p1, *p2, *p3, *p1a;

  r = 0;

  for (i = 0; i < m_CgxDataSize; i++) {
    /* Continuity check of the same data */
    for (j = 1; (i + j) < m_CgxDataSize; j++) {
      if (m_aCgxData[i] != m_aCgxData[i + j]) {
        break;
      }
    }

    m_algorithm[i] = RUNLENGTH; /* Compression Algorithms */
    m_compLength[i] = j;        /* Data length to be compressed */

    /* Investigation of alternating continuity of two data */
    if (((i + 3) < m_CgxDataSize)
        && (m_aCgxData[i] == m_aCgxData[i + 2])
        && (m_aCgxData[i + 1] == m_aCgxData[i + 3])) {
      for (j = 2; ((i + j) < (m_CgxDataSize - 1)); j += 2) {
        if ((m_aCgxData[i + 0] != m_aCgxData[i + j + 0]) ||
            (m_aCgxData[i + 1] != m_aCgxData[i + j + 1])) {
          break;
        }
      }

      if (m_aCgxData[i] == m_aCgxData[i + j]) {
        j++;
      }
    }

    m_algorithm[i] = RUNLENGTH2; /* Compression Algorithms */
    m_compLength[i] = j;         /* Data length to be compressed */

    /* 0 continuity investigation */
    if (m_aCgxData[i] == 0) {
      for (j = 1; (i + j) < m_CgxDataSize; j++) {
        if (m_aCgxData[i + j] != 0) {
          break;
        }
      }

      m_algorithm[i] = INCRUNLEN; /* Compression Algorithms */
      m_compLength[i] = j;        /* Data length to be compressed */
    }

    if (m_compLength[i] >= LONG_ENOUGH) {
      i += (m_compLength[i] - 1);
      continue;
    }

    /* Investigating references to previously occurring data */
    l = 0;
    hn = hashf(m_aCgxData[i + 0], m_aCgxData[i + 1],
               m_aCgxData[i + 2], m_aCgxData[i + 3]);

    while ((j = scanHash(hn, NREF)) != EMPTY) {
      p1 = p1a = &m_aCgxData[i];
      p2 = &m_aCgxData[j];
      p3 = m_aCgxData + m_CgxDataSize;

      while ((*p1 == *p2) && (p1 < p3)) {
        ++p1, ++p2;
      }

      if ((k = p1 - p1a) > l) {
        l = k, r = j; /* The DAT scanned currently has a longer match */
      } else if ((k == l) && (j > r)) {
        l = k, r = j; /* If identical lengths, proxy takes precedence */
      }
    }

    if (l >= 3 && m_compLength[i] < l) {
      m_algorithm[i] = REFERENCE; /* Compression Algorithms */
      m_compLength[i] = l;        /* Data length to be compressed */
      m_reference[i] = r;         /* Reference Position */
    }

    /* Check for L-R flipped references to previously occurring data */
    l = 0;
    hh = hashf(m_invertTable[m_aCgxData[i + 0]],
               m_invertTable[m_aCgxData[i + 1]],
               m_invertTable[m_aCgxData[i + 2]],
               m_invertTable[m_aCgxData[i + 3]]);

    while ((j = scanHash(hh, NREF)) != EMPTY) {
      p1 = p1a = m_aCgxData + i;
      p2 = m_aCgxData + j;
      p3 = m_aCgxData + m_CgxDataSize;

      while ((*p1 == m_invertTable[*p2]) && (p1 < p3)) {
        ++p1, ++p2;
      }

      if ((k = p1 - p1a) > l) {
        l = k, r = j; /* The data scanned ATM has a longer match */
      } else if ((k == l) && (j > r)) {
        l = k, r = j; /* If identical lengths, proxy takes precedence */
      }
    }

    if (l >= 3 && m_compLength[i] < l) {
      m_algorithm[i] = REFERLR; /* Compression Algorithms */
      m_compLength[i] = l;      /* Data length to be compressed */
      m_reference[i] = r;       /* Reference Position */
    }

    /* Investigating U/D references to previously occurring data */
    l = 0;

    while ((j = scanHash(hn, VREF)) != EMPTY) {
      p1 = p1a = &m_aCgxData[i];
      p2 = &m_aCgxData[j];
      p3 = m_aCgxData + m_CgxDataSize;

      while ((*p1 == *p2) && (p2 >= m_aCgxData) && (p1 < p3)) {
        ++p1, --p2;
      }

      if ((k = p1 - p1a) > l) {
        l = k, r = j; /* The data scanned ATM has a longer match */
      } else if ((k == l) && (j > r)) {
        l = k, r = j; /* If identical lengths, proxy takes precedence */
      }
    }

    if ((l >= 3) && (m_compLength[i] < l)) {
      m_algorithm[i] = REFERUD; /* Compression Algorithms */
      m_compLength[i] = l;      /* Data length to be compressed */
      m_reference[i] = r;       /* Reference Position */
    }

    registNHash(i, hn);
    registVHash(i);

    if (m_compLength[i] >= LONG_ENOUGH) {
      i += (m_compLength[i] - 1);
      continue;
    }
  }

  return;
}

static u8 reverseOctet(u8 byte)
{
  u8 tmp;
  int i;

  for (i = 0, tmp = 0; i < 8; i++) {
    tmp <<= 1;
    tmp |= ((byte >> i) & 1);
  }

  return tmp;
}

static void init_pack(void)
{
  u8 *p = &m_invertTable[0], b = 0;

  while (p < &m_invertTable[256]) { *p++ = reverseOctet(b++); }

  memset(m_nrefHash, EMPTY, HASHTABLESIZE - 1);
  memset(m_vrefHash, EMPTY, HASHTABLESIZE - 1);
  return;
}

static void pack(void)
{
  /* Initial settings required for compression */
  m_packedSize = 0;
  init_pack();
  checkData();
  judge();
  outData();
  return;
}

#define D_RUNLENGTH  0x20 /*  32 : 001 */
#define D_RUNLENGTH2 0x40 /*  64 : 010 */
#define D_INCRUNLEN  0x60 /*  96 : 011 */
#define D_REFERENCE  0x80 /* 128 : 100 */
#define D_REFERLR    0xA0 /* 160 : 101 */
#define D_REFERUD    0xC0 /* 192 : 110 */
#define D_SPECIAL    0xE0 /* 224 : 111 */
static void unpack(void)
{    
  u8 *src, *bin_src, *dest, DataType;
  u16 len;
  int CheckDataCount;

  src = bin_src = (u8 *)malloc(sizeof(u8) * v_size(&m_aBinData));
  dest = &m_aCgxData[0];
  CheckDataCount = 0;
  memcpy(src, m_aBinData.d, v_size(&m_aBinData));
  memset(m_aCgxData, 0, CGB_VRAM_SIZE);

  /* Extracting Header Information */
  while (*src != ENDMARK) {
    if ((*src & D_SPECIAL) == D_SPECIAL) {
      /* 2-byte format */
      DataType = (*src << 3) & D_SPECIAL;
      len = ((*src & 3) << 8) + *(src + 1) + 1;
      src += 2;
    }
    else {
      /* 1-byte format */
      DataType = *src & D_SPECIAL;
      len = (u16)(*src & 31) + 1;
      src++;
    }
        
    if ((DataType & D_REFERENCE) != 0) {
      /* It is reference data */
      u8 *src2 = src;

      if ((*src & D_REFERENCE) != 0) {
        src2 = dest - ((*src & 127) + 1);
        src++;
      }
      else {
        src2 = &m_aCgxData[0] + (*src << 8) + *(src + 1);
        src+=2;
      }

      if (DataType == D_REFERENCE) {
        /* Regular reference expansion */
        while (len-- > 0) {
          *dest++ = *src2++;
          CheckDataCount++;
        }
      }
      else if (DataType == D_REFERLR) {
        /* Left-right reversed reference expansion */
        while (len-- > 0) {
          *dest++ = reverseOctet(*src2++);
          CheckDataCount++;
        }
      }
      else if (DataType == D_REFERUD) {
        /* Upside-down reference expansion */
        while (len-- > 0) {
          *dest++ = *src2--;
          CheckDataCount++;
        }
      }
    }
    else if (DataType == D_RUNLENGTH) {
      /* Expanding a series of data */
      while (len-- > 0) {
        *dest++ = *src;
        CheckDataCount++;
      }
      src++;
    }
    else if (DataType == D_RUNLENGTH2) {
      /* Two data alternating successive developments */
      while (len-- > 0) {
        *dest++ = src[0];
        CheckDataCount++;

        if (len-- > 0) {
          *dest++ = src[1];
          CheckDataCount++;
        }
        else {
          break;
        }
      }

      src += 2;
    }
    else if (DataType == D_INCRUNLEN) {
      /* +1 Data sequence expansion */
      while (len-- > 0) {
        *dest++ = 0;
        CheckDataCount++;
      }
    }
    else {
      /* Decompressing uncompressed data */
      while (len-- > 0) {
        *dest++ = *src++;
        CheckDataCount++;
      }
    }
  }

  m_CgxDataSize = CheckDataCount;
  free(bin_src);
  bin_src = NULL;
  return;
}

static void usage(void)
{
  printf("\nUsage:\n\tvpack [mode] <in_file> <out_file>\n");
  printf("\nMode:\n\tP : Pack\n\tU : Unpack\n");
  return;
}

int main(int argc, char *argv[])
{
  FILE *ifile = NULL, *ofile = NULL;
  char *s;

  if (argc != 4) {
    usage();
    exit(EXIT_FAILURE);
  }

  if ((s = argv[1], s[1] || strpbrk(s, "PUpu") == NULL) ||
      (s = argv[2], (ifile = fopen(s, "rb")) == NULL) ||
      (s = argv[3], (ofile = fopen(s, "wb")) == NULL)) {
    printf("??? %s\n", s);
    usage();
    exit(EXIT_FAILURE);
  }

  fseek(ifile, 0, SEEK_END);
  m_CgxDataSize = ftell(ifile);
  rewind(ifile);

  if (m_CgxDataSize > CGB_VRAM_SIZE) {
    printf("[VAL] %u > [MAX] %u\n", m_CgxDataSize, CGB_VRAM_SIZE);
    fclose(ofile);
    unlink(argv[3]);
    fclose(ifile);
    exit(EXIT_FAILURE);
  }

  if (toupper(*argv[1]) == 'P') {
    fread(m_aCgxData, sizeof(u8), m_CgxDataSize, ifile);
    v_init(&m_aBinData, sizeof(u8));
    pack();
    fwrite(m_aBinData.d, sizeof(u8), m_aBinData.s*m_aBinData.es, ofile);
    v_free(&m_aBinData);    
  }
  else {
    int c;

    v_init(&m_aBinData, sizeof(u8));

    while ((c = fgetc(ifile)) != EOF) { v_push(&m_aBinData, &c); }

    unpack();
    fwrite(&m_aCgxData[0], sizeof(u8), m_CgxDataSize, ofile);
    v_free(&m_aBinData);
  }

  fflush(ofile);
  fclose(ofile);
  fclose(ifile);
  exit(EXIT_SUCCESS);
}





/* ------------------------------------------------------------------ */
/*                                                                    */
/*                   HISTORICAL PRESERVATION OBJECT                   */
/*                                                                    */
/* ------------------------------------------------------------------ */
