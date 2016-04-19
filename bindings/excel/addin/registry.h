#pragma once

namespace prxl { 
    
    /// Contains environment-related functions
    namespace env {

#undef LONG

    /// Reads from system registry a Premia Excel setting.
    /// \param branch - a branch to be read
    /// \param key - a key in the branch to be read
    /// \return the registry key value
    inline _bstr_t getSetting(_bstr_t const & branch, _bstr_t const & key)
    {
        HKEY hKey;
        wchar_t value[10000];
        DWORD dwBufLen=10000;
        long lRet;

        lRet = RegOpenKeyEx( HKEY_CURRENT_USER,
            _bstr_t("Software\\VB and VBA Program Settings\\Premia_9\\") + branch,
            0, KEY_QUERY_VALUE, &hKey );
        if( lRet != ERROR_SUCCESS )
            return "-=<Error: Registry branch is not found>=-";

        lRet = RegQueryValueEx( hKey, key, NULL, NULL,
            (LPBYTE) value, &dwBufLen);
        if( (lRet != ERROR_SUCCESS) || (dwBufLen > 10000) )
            return "-=<Error: Registry key is not found>=-";

        RegCloseKey( hKey );

        return value;
    }

    /// \return path to root of Premia helpfiles
    inline _bstr_t getManPath()
    {
        static _bstr_t x = getSetting("Man", "Path");
        return x;
    }

    inline _bstr_t getDataPath()
    {
        static _bstr_t x = getSetting("data", "path");
        return x;
    }

    /// \return path to a directory where to locate a XLL with Premia add-in formulas
    inline _bstr_t getXLLPath()
    {
        static _bstr_t x = getSetting("XLL", "Path");
        return x;
    }

    /// \return path to a directory where to locate a XLA with old Premia 8 interface
    inline _bstr_t getXLAPath()
    {
        static _bstr_t x = getSetting("XLA", "Path");
        return x;
    }

    inline _bstr_t currentDirectory()
    {
        wchar_t buffer[1000];
        ::GetCurrentDirectory(1000, buffer);
        return buffer;
    }

} }