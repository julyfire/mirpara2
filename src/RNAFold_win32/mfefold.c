#include <jni.h>
#include "jmipara2_MfeFold.h"
#include  <math.h>
#include  <string.h>
#include  "utils.h"
#include  "fold_vars.h"
#include  "fold.h"


jfieldID structure;

JNIEXPORT void JNICALL Java_jmipara2_MfeFold_initIDs
  (JNIEnv *env, jclass cls)
{
	structure=(*env)->GetStaticFieldID(env,cls,"structure","Ljava/lang/String;");
    if(structure==NULL) return;
}

JNIEXPORT jfloat JNICALL Java_jmipara2_MfeFold_fold
  (JNIEnv *env, jclass cls, jstring ss, jfloat t)
{
	jstring str;
	float e;
	char* seq;
	char* struc;
	
	seq=(char*)(*env)->GetStringUTFChars(env,ss,NULL);
		if(seq==NULL) return 0;
	/*char* struc=(char*)space(sizeof(char)*(strlen(seq)+1));*/
    str=(*env)->GetStaticObjectField(env,cls,structure);
    struc=(char*)(*env)->GetStringUTFChars(env,str,NULL);
                /*if(struc==NULL) return;*/
	temperature = t;
    do_backtrack=0;

	/*initialize_fold(strlen(seq));*/
	update_fold_params();
        
    /*if the structure has been provided, calculate the energy of the structure*/
    if(strlen(struc)){
        dangles=0;
        update_fold_params();
        e=energy_of_struct(seq,struc);
    }
    /*or to fold the sequence to get the structure and energy*/
    else{
        struc=(char*)space(sizeof(char)*(strlen(seq)+1));
        e=fold(seq,struc);
        str=(*env)->NewStringUTF(env,struc);
		if(str==NULL) return 0;
        (*env)->SetStaticObjectField(env,cls,structure,str);/*return structure*/
        free_arrays();/*frees the memory allocated by fold()*/
    }

	
	(*env)->ReleaseStringUTFChars(env,ss,seq);
    /*(*env)->ReleaseStringUTFChars(env,str,struc);*/

	return e;

}
