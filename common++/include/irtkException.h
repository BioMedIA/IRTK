/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKEXCEPTION_H

#define _IRTKEXCEPTION_H

#include <string>

#include <irtkObject.h>

/** Exception class for the irtk library. */
class irtkException : public irtkObject
{
public:
  /** Constructor takes message, file name and line number at which exception
      was thrown.
      \param message The message.
      \param file The file name. You can use the __FILE__ macro to get the
      file name.
      \param line The line number at which the exception was thrown. You can
      use the __LINE__ macro to get the line number. */
  irtkException(const std::string& message, const std::string& fileName = "",
                unsigned int line = 0);

  /** Returns the message. */
  const std::string& GetMessage() const;

  /** Returns the name of the file in which the exception was thrown. */
  const std::string& GetFileName() const;

  /** Returns the line number at which the exception was thrown. */
  unsigned int GetLine() const;

  /** Outputs the exception to a stream.
      \param os The output stream.
      \param ex The exception. */
  friend std::ostream& operator<<(std::ostream& os, const irtkException& ex);

protected:
  /** The message. */
  std::string _Message;

  /** The file name. */
  std::string _FileName;

  /** The line number at which the exception was thrown. */
  unsigned int _Line;
};

inline irtkException::irtkException(const std::string& strMessage,
                                    const std::string& fileName, unsigned int line)
{
  _Message = strMessage;
  _FileName = fileName;
  _Line = line;
}

inline const std::string& irtkException::GetMessage() const
{
  return _Message;
}

inline const std::string& irtkException::GetFileName() const
{
  return _FileName;
}

inline unsigned int irtkException::GetLine() const
{
  return _Line;
}

#endif // _IRTKEXCEPTION_H

