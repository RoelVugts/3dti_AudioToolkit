/**
* \class CBuffer
*
* \brief Declaration of CBuffer interface.
* \date	July 2016
*
* \authors 3DI-DIANA Research Group (University of Malaga), in alphabetical order: M. Cuevas-Rodriguez, C. Garre,  D. Gonzalez-Toledo, E.J. de la Rubia-Cuestas, L. Molina-Tanco ||
* Coordinated by , A. Reyes-Lecuona (University of Malaga) and L.Picinali (Imperial College London) ||
* \b Contact: areyes@uma.es and l.picinali@imperial.ac.uk
*
* \b Contributions: (additional authors/contributors can be added here)
*
* \b Project: 3DTI (3D-games for TUNing and lEarnINg about hearing aids) ||
* \b Website: http://3d-tune-in.eu/
*
* \b Copyright: University of Malaga and Imperial College London - 2018
*
* \b Licence: This copy of 3dti_AudioToolkit is licensed to you under the terms described in the 3DTI_AUDIOTOOLKIT_LICENSE file included in this distribution.
*
* \b Acknowledgement: This project has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement No 644051
*/

#ifndef _CBUFFER_H_
#define _CBUFFER_H_

#define _USE_MATH_DEFINES // TODO: Test in windows! Might also be problematic for other platforms??
#include <cmath>
#include <vector>
#include <algorithm>
// #include <cassert>
#include <Common/ErrorHandler.h>
// #include <Common/Magnitudes.h>
//#include <initializer_list>

/*! \file */

namespace Common {

	/** \details This is a template class to manage audio streamers and buffers
	*/
	template <
		unsigned int NChannels,
		class stored
	>
		class CBuffer : public std::vector<stored, std::allocator<stored>>
	{
	public:
		using std::vector<stored>::vector;    //   inherit all std::vector constructors
		using std::vector<stored>::size;      //   MacOSX clang seems to need this to compile
		using std::vector<stored>::begin;     //   MacOSX clang seems to need this to compile
		using std::vector<stored>::end;       //   MacOSX clang seems to need this to compile
		using std::vector<stored>::resize;     //   MacOSX clang seems to need this to compile

		/** \brief Get number of channels in the buffer
		*	\retval nChannels number of channels
		*   \eh Nothing is reported to the error handler.
		*/
		constexpr unsigned int GetNChannels() const
		{
			return NChannels;
		}

		/** \brief Add all sample values from another buffer
		*   \eh Nothing is reported to the error handler.
		*/
		CBuffer<NChannels, stored> & operator+= (const CBuffer<NChannels, stored> & oth)
		{
			//assert(GetNChannels()==oth.GetNChannels()); // TODO: agree on error handling
			//assert(size()==oth.size());
			//ASSERT(GetNChannels() == oth.GetNChannels(), "Attempt to mix two buffers of different sizes", ""); 
			//ASSERT(size() == oth.size(), "Attempt to mix two buffers of different sizes", "");

			std::transform(begin(), end(), oth.begin(), begin(), [](stored a, stored b) { return a + b; });
			return *this;
		}

		/** \brief Substract all sample values of another buffer
		*   \eh Nothing is reported to the error handler.
		*/
		CBuffer<NChannels, stored> & operator-= (const CBuffer<NChannels, stored> & oth)
		{
			//assert(GetNChannels()==oth.GetNChannels()); // TODO: agree on error handling
			//assert(size()==oth.size());
			//ASSERT(GetNChannels() == oth.GetNChannels(), "Attempt to mix two buffers of different sizes", ""); 
			//ASSERT(size() == oth.size(), "Attempt to mix two buffers of different sizes", "");

			std::transform(begin(), end(), oth.begin(), begin(), [](stored a, stored b) { return a - b; });
			return *this;
		}

		/** \brief Multiply the values in the buffer by a constant
		*   \eh Nothing is reported to the error handler.
		*/
		// Following recommendations for binary arithmetic operators in: http://en.cppreference.com/w/cpp/language/operators
		friend CBuffer<NChannels, stored> operator* (CBuffer<NChannels, stored> buffer, const stored& gain)
		{
			buffer.ApplyGain(gain);
			return buffer;
		}

		/** \brief Add the values of two buffers
		*   \eh Nothing is reported to the error handler.
		*/
		// Following recommendations for binary arithmetic operators in: http://en.cppreference.com/w/cpp/language/operators
		friend CBuffer<NChannels, stored> operator+ (CBuffer<NChannels, stored> lhs, const CBuffer<NChannels, stored>& rhs)
		{
			lhs += rhs;
			return lhs;
		}

		/** \brief Multiply the values in the buffer by a constant gain
		*	\param [in] gain constant gain value to apply
		*   \eh Nothing is reported to the error handler.
		*/
		void ApplyGain(stored gain)
		{
			//SET_RESULT(RESULT_OK, "Gain applied to buffer succesfully");
			std::for_each(begin(), end(), [gain](stored & a) { a *= gain; });
		}

		/** \brief Fill nToFill samples of the buffer with the same value
		*	\param [in] nToFill number of samples to fill
		*	\param [in] value value to fill with
		*   \eh On success, RESULT_OK is reported to the error handler.
		*/
		void Fill(size_t nToFill, stored value)
		{
			SET_RESULT(RESULT_OK, "Buffer filled with single value succesfully");
			// Assign is the fastest implementation, after memset. See: http://stackoverflow.com/questions/8848575/fastest-way-to-reset-every-value-of-stdvectorint-to-0
			this->assign(nToFill, value);
		}

		/** \brief Get number of samples in each channel of the buffer
		*	\pre Buffer must have at least one channel
		*	\retval nSamples number of samples per channel
		*   \eh Nothing is reported to the error handler.
		*/
		unsigned long GetNsamples() const
		{
			return size() / GetNChannels();
		}

		/** \brief Interlace two mono buffers into one stereo buffer
		*	\param [in] left mono buffer for left channel
		*	\param [in] right mono buffer for right channel
		*	\pre this must be a stereo buffer
		*	\pre left and right must have the same size
		*   \eh On error, an error code is reported to the error handler.
		*/
		void Interlace(CBuffer<1, stored> & left, CBuffer<1, stored> & right)
		{
			// Preconditions check
			ASSERT(GetNChannels() == 2, RESULT_ERROR_BADSIZE, "Attempt to interlace into a non-stereo buffer", "");
			ASSERT(left.size() == right.size(), RESULT_ERROR_BADSIZE, "Attempt to interlace two mono buffers of different length", "");
			//SET_RESULT(RESULT_OK, "Stereo buffer interlaced from two mono buffers succesfully");

			// Start with a clean buffer
			this->clear();
			this->reserve(2 * left.size());

			// Interlace channels
			for (int sample = 0; sample < left.size(); sample++)
			{
				this->push_back(left[sample]);
				this->push_back(right[sample]);
			}
		}

		/** \brief Set buffer from a pure tone
		*	\details For debugging purposes
		*	\param [in] samplingRate sampling rate in Hzs
		*	\param [in] frequency frequency of tone in Hzs
		*	\param [in] amplitude amplitude of tone, as gain (typically, between 0.0 and 1.0)
		*	\param [in] phase phase of tone in radians
		*   \eh Nothing is reported to the error handler.
		*/
		void SetFromTone(float samplingRate, float frequency, float amplitude = 1.0f, float phase = 0.0f)
		{
			for (size_t i = 0; i < GetNsamples(); i++)
			{
				for (int c = 0; c < GetNChannels(); c++)
				{
					float value = amplitude * std::sin(2.0f * M_PI*frequency*((float)i / samplingRate) + phase);
					(*this)[i] = value;
					//WATCH(WV_BUFFER_TEST, value, float);
				}
			}
		}

		/** \brief Calculate power of a mono buffer signal
		*	\retval power power
		*   \eh Nothing is reported to the error handler.
		*/
		float GetPower()
		{
			return GetAutocorrelation( 0 );
		}

		/** \brief Calculate autocorrelation of a mono buffer 
		*	\param [in] n displacement in number of samples 
		*	\details The output is not normalized. To normalize, divide it by GetPower()
		*	\retval coefficient n-th coefficient of autocorrelation
		*   \eh Nothing is reported to the error handler.
		*/
		float GetAutocorrelation( int n )
		{
			ASSERT(GetNChannels() == 1, RESULT_ERROR_BADSIZE, "Attempt to calculate autocorrelation of a non-mono buffer", "");
			ASSERT((*this).size() > 0, RESULT_ERROR_BADSIZE, "Attempt to calculate autocorrelation of a empty buffer", "");
			ASSERT((*this).size() > n && n >= 0, RESULT_ERROR_INVALID_PARAM, "Invalid displacement in GetAutocorrelation", "");

			int bufferSize = (*this).size();
			float acum = 0.0f;

			int overlapedSamples = bufferSize - n;

			// Calculate Autocorrelation = (1/(bufferSize-n)) *sum(x[i] * x[i-n])
			for(int c=0; c < overlapedSamples; c++ )
				acum += (*this)[c] * (*this)[c+n];
			
			return overlapedSamples <= 0 ? 0 : acum / ((float)overlapedSamples);
		}
	};
}

/** \brief One channel specialization of CBuffer
*/
template<class stored>
using CMonoBuffer = Common::CBuffer<1,stored>;

/** \brief Two channels specialization of CBuffer
*/
template<class stored>
using CStereoBuffer = Common::CBuffer<2,stored>;

/** \brief Non-enforcing buffer 
*	\details Current implementation does not inherits from CBuffer
*/
template<class stored>
using CMultiChannelBuffer = std::vector<stored>; // TO DO: why isnt it CBuffer?????

#endif
