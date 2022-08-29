/*
Copyright (c) 2022, Jeff Hurchalla

Original source file prior to modifications was:
https://github.com/bbuhrow/yafu/blob/25b65990d6501b0a71e69963fb59c1fc4ab28df1/factor/gmp-ecm/microecm.c

IMPORTANT
Currently (and temporarily) this file and all of its modifications from the
original file are unavailable under any license.  For now, you are explicitly
*NOT* allowed to use, copy, distribute, modify, or share this software for any
purpose, without permission from the author.

You can expect this file will soon have a normal permissive license.
*/

/*
The following is reproduced to comply with the original source file:

Copyright (c) 2014, Ben Buhrow
All rights reserved.

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
*/



#ifndef MICROECM_GETFACTOR_ECM_C_API_H_INCLUDED
#define MICROECM_GETFACTOR_ECM_C_API_H_INCLUDED

#ifdef __cplusplus   // C compilers skip this ifdef section
extern "C" {
#endif

// This is the C API for ECM factoring.
// getfactor_ecm() returns a *single* factor.

// getfactor_ecm() returns 1 if unable to find a factor of n,
// otherwise returns a factor of n.
//
// Prior to your first call of getfactor_ecm(), set *pran = 0 (or some arbitrary
// value); after that, don't change *pran.
// FYI: *pran is used within microecm_c.c by a random number generator, and
// holds the current value of a pseudo random sequence.  Your first assigment
// to *pran seeds the sequence, and after seeding it you don't want to
// change *pran, since that would restart the sequence.

uint64_t getfactor_ecm(uint64_t n, uint64_t *pran);


#ifdef __cplusplus
}
#endif

#endif
