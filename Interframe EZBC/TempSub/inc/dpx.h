#include "basic.h"
#define Gamma 2.2
#define Refblack  95            //95
#define Refwhite  685           //685
//#define KODAK
#define PEAK 1023
#define DCtest

/*generic image data header */
typedef struct file_information
{
  U32 magic_num;                /* magic number 0x53445058 (SDPX) or 0x58504453 (XPDS) */
  U32 offset;                   /* offset to image data in bytes */
  ASCII vers[8];                /* which header format version is being used (v1.0) */
  U32 file_size;                /* file size in bytes */
  ASCII Undefined1[16];
  ASCII Undefined2[624];
  ASCII Undefined3[108];

/*U32   ditto_key; *//* read time short cut - 0 = same, 1 = new */
/*U32   gen_hdr_size; *//* generic header length in bytes */
/*U32   ind_hdr_size; *//* industry header length in bytes */
/*U32   user_data_size; *//* user-defined data length in bytes */
/*ASCII file_name[100]; *//* iamge file name */
/*ASCII create_time[24]; *//* file creation date "yyyy:mm:dd:hh:mm:ss:LTZ" */
/*ASCII creator[100]; *//* file creator's name */
/*ASCII project[200]; *//* project name */
/*ASCII copyright[200]; *//* right to use or copyright info */
/*U32   key; *//* encryption ( FFFFFFFF = unencrypted ) */
/*ASCII Reserved[104]; *//* reserved field TBD (need to pad) */
}
FileInformation;


/*Image Information Header*/
typedef struct _image_information
{
  U16 orientation;              /* image orientation */
  U16 element_number;           /* number of image elements */
  U32 pixels_per_line;          /* or x value */
  U32 lines_per_image_ele;      /* or y value, per element */
  ASCII Undefined1[16];

  /*struct _image_element */
  /*{ */
  U32 data_sign;                /* data sign (0 = unsigned, 1 = signed ) */

  /* "Core set images are unsigned" */
/*U32    ref_low_data; *//* reference low data code value */
/*R32    ref_low_quantity; *//* reference low quantity represented */
/*U32    ref_high_data; *//* reference high data code value */
/*R32    ref_high_quantity; *//* reference high quantity represented */
  U8 descriptor;                /* descriptor for image element */
  U8 transfer;                  /* transfer characteristics for element */
  U8 colorimetric;              /* colormetric specification for element */
  U8 bit_size;                  /* bit size for element */
  U16 packing;                  /* packing for element */
  U16 encoding;                 /* encoding for element */
  U32 data_offset;              /* offset to data of element */


/*U32    eol_padding; *//* end of line padding used in element */
/*U32    eo_image_padding; *//* end of image padding used in element */
/*ASCII  description[32]; *//* description of element */
/*} image_element[8]; *//* NOTE THERE ARE EIGHT OF THESE */
#ifdef DCtest
  ASCII Undefined3[7380];
#else
  ASCII Undefined3[64724];
#endif
/*U8 reserved[52]; *//* reserved for future use (padding) */
}
Image_Information;

/*Image Orientation Information*/
typedef struct _image_orientation
{
  U32 x_offset;                 /* X offset */
  U32 y_offset;                 /* Y offset */
  R32 x_center;                 /* X center */
  R32 y_center;                 /* Y center */
  U32 x_orig_size;              /* X original size */
  U32 y_orig_size;              /* Y original size */
  ASCII file_name[250];         /* source image file name */
  ASCII creation_time[24];      /* source image creation date and time */
  ASCII input_dev[32];          /* input device name */
  ASCII input_serial[32];       /* input device serial number */
  U16 border[4];                /* border validity (XL, XR, YT, YB) */
  U32 pixel_aspect[2];          /* pixel aspect ratio (H:V) */
  U8 reserved[28];              /* reserved for future use (padding) */
}
Image_Orientation;


typedef struct dpx
{
  FileInformation fileinfo;
  Image_Information imageinfo;
}
DPXHEADER;


/*Image Data
The image data is stored as an array of 32-bit elements made up of four signed or unsigned character values. Unsigned values are the default for DPX image data. */
typedef struct _image_data_element
{
  U32 *Data;
}
Image_Data_Element;




U16 *read_dpx( char *dpxname, struct dpx **dpxheader, int *nr, int *nc );
void write_dpx( char *dpxname, float *RGB, videoinfo info );
void dpxrange( U16 * RGBframe, int row, int col, U16 * Bmin, U16 * Bmax,
               U16 * Gmin, U16 * Gmax, U16 * Rmin, U16 * Rmax );

float *linearize( U16 * RGBframe, int row, int col, videoinfo info );
float *Nonlinearize( float *LRGB, int row, int col, videoinfo info );
void Pack( U32 * data, int comp, U8 bitdepth );
void save_dpx_header( videoinfo info, int index );
struct dpx *read_dpx_header( char *filename );
void compute_colorLUT( int bits, enum CODING_DOMAIN coding_domain );
float convertLOG( int peak, float data, enum CODING_DOMAIN coding_domain );
