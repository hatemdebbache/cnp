Attribute VB_Name = "SYSLIN"
Option Explicit
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << SYSLIN >> version 03
'  Requis   :
'  Externes :
'
'  Descrip. : Résolution directe de systèmes linéaires
'             du type Ax = b
'
' ================================================
'
'            Système à matrice tridiagonale
'
' ================================================
Function Thomas(ByVal n As Integer, _
                ByRef d() As Double, _
                ByRef L() As Double, _
                ByRef U() As Double, _
                ByRef b() As Double, _
                ByRef x() As Double) As Integer

    Dim k As Integer
    Dim c() As Double, p As Double
    
    ReDim c(1 To n), x(1 To n)
    
    If n < 2 Then Stop ' Données non valides
    p = d(1)
    If (p = 0#) Then Thomas = 0: Exit Function ' Singularité
    x(1) = b(1) / p
    For k = 2 To n
        c(k) = U(k - 1) / p
        p = d(k) - L(k) * c(k)
        If (p = 0) Then Thomas = 0: Exit Function
        x(k) = (b(k) - L(k) * x(k - 1)) / p
    Next k
    For k = n - 1 To 1 Step -1
        x(k) = x(k) - c(k + 1) * x(k + 1)
    Next k
    Thomas = 1 ' OK
End Function
Function Thomas2(ByVal n As Integer, _
                 ByRef d() As Double, _
                 ByRef L() As Double, _
                 ByRef U() As Double, _
                 ByRef b() As Double) As Integer

    Dim k As Integer
    Dim c() As Double, p As Double
    
    ReDim c(1 To n)
    
    If n < 2 Then Stop ' Données non valides
    p = d(1)
    If (p = 0#) Then Thomas2 = 0: Exit Function ' Singularité
    b(1) = b(1) / p
    For k = 2 To n
        c(k) = U(k - 1) / p
        p = d(k) - L(k) * c(k)
        If (p = 0) Then Thomas2 = 0: Exit Function
        b(k) = (b(k) - L(k) * b(k - 1)) / p
    Next k
    For k = n - 1 To 1 Step -1
        b(k) = b(k) - c(k + 1) * b(k + 1)
    Next k
    Thomas2 = 1 ' OK
End Function
Function ThomasDiag(ByVal n As Integer, _
                    ByVal d As Double, _
                    ByRef L() As Double, _
                    ByRef U() As Double, _
                    ByRef b() As Double, _
                    ByRef x() As Double) As Integer

    Dim k As Integer
    Dim c() As Double, p As Double
    
    ReDim c(1 To n), x(1 To n)
    
    If n < 2 Then Stop ' Données non valides
    p = d
    If (p = 0#) Then ThomasDiag = 0: Exit Function ' Singularité
    x(1) = b(1) / p
    For k = 2 To n
        c(k) = U(k - 1) / p
        p = d - L(k) * c(k)
        If (p = 0) Then ThomasDiag = 0: Exit Function
        x(k) = (b(k) - L(k) * x(k - 1)) / p
    Next k
    For k = n - 1 To 1 Step -1
        x(k) = x(k) - c(k + 1) * x(k + 1)
    Next k
    ThomasDiag = 1 ' OK
End Function
Function ThomasDiag2(ByVal n As Integer, _
                     ByVal d As Double, _
                     ByRef L() As Double, _
                     ByRef U() As Double, _
                     ByRef b() As Double) As Integer

    Dim k As Integer
    Dim c() As Double, p As Double
    
    ReDim c(1 To n)
    
    If n < 2 Then Stop ' Données non valides
    p = d
    If (p = 0#) Then ThomasDiag2 = 0: Exit Function ' Singularité
    b(1) = b(1) / p
    For k = 2 To n
        c(k) = U(k - 1) / p
        p = d - L(k) * c(k)
        If (p = 0) Then ThomasDiag2 = 0: Exit Function
        b(k) = (b(k) - L(k) * x(k - 1)) / p
    Next k
    For k = n - 1 To 1 Step -1
        b(k) = b(k) - c(k + 1) * b(k + 1)
    Next k
    ThomasDiag2 = 1 ' OK
End Function
Function ThomasConst(ByVal n As Integer, _
                     ByVal d As Double, _
                     ByVal L As Double, _
                     ByVal U As Double, _
                     ByRef b() As Double, _
                     ByRef x() As Double) As Integer

    Dim k As Integer
    Dim c() As Double, p As Double
    
    ReDim c(1 To n), x(1 To n)
    
    If n < 2 Then Stop ' Données non valides
    p = d
    If (p = 0#) Then ThomasConst = 0: Exit Function ' Singularité
    x(1) = b(1) / p
    For k = 2 To n
        c(k) = U / p
        p = d - L * c(k)
        If (p = 0) Then ThomasConst = 0: Exit Function
        x(k) = (b(k) - L * x(k - 1)) / p
    Next k
    For k = n - 1 To 1 Step -1
        x(k) = x(k) - c(k + 1) * x(k + 1)
    Next k
    ThomasConst = 1 ' OK
End Function
Function ThomasConst2(ByVal n As Integer, _
                      ByVal d As Double, _
                      ByVal L As Double, _
                      ByVal U As Double, _
                      ByRef b() As Double) As Integer

    Dim k As Integer
    Dim c() As Double, p As Double
    
    ReDim c(1 To n)
    
    If n < 2 Then Stop ' Données non valides
    p = d
    If (p = 0#) Then ThomasConst2 = 0: Exit Function ' Singularité
    b(1) = b(1) / p
    For k = 2 To n
        c(k) = U / p
        p = d - L * c(k)
        If (p = 0) Then ThomasConst2 = 0: Exit Function
        b(k) = (b(k) - L * b(k - 1)) / p
    Next k
    For k = n - 1 To 1 Step -1
        b(k) = b(k) - c(k + 1) * b(k + 1)
    Next k
    ThomasConst2 = 1 ' OK
End Function
' ================================================
'
'            Système à matrice pleine
'
' ================================================
Function LU_FacSol(ByVal n As Integer, _
                   ByRef a() As Double, _
                   ByRef b() As Double, _
                   ByRef x() As Double, _
          Optional ByVal Cdet As Boolean = False) As Double
   
    Dim s As Integer, det As Double
    Dim ndx() As Integer
   
    det = LU_Fac(n, a(), ndx(), s)
    If det = 1 Then
        Call LU_Sol(n, a(), ndx(), b(), x())
        If Cdet Then det = LU_Det(n, a(), ndx(), s)
    End If
    LU_FacSol = det
End Function
Function LU_FacSol2(ByVal n As Integer, _
                    ByRef a() As Double, _
                    ByRef b() As Double, _
           Optional ByVal Cdet As Boolean = False) As Double
   
    Dim s As Integer, det As Double
    Dim ndx() As Integer
   
    det = LU_Fac(n, a(), ndx())
    If det = 1 Then
        Call LU_Sol2(n, a(), ndx(), b())
        If Cdet Then det = LU_Det(n, a(), ndx(), det)
    End If
    LU_FacSol2 = det
End Function
Function LU_FacInv(ByVal n As Integer, _
                   ByRef a() As Double, _
                   ByRef INV() As Double, _
          Optional ByVal Cdet As Boolean = False) As Double

    Dim s As Integer, det As Double
    Dim ndx() As Integer
    
    det = LU_Fac(n, a(), ndx())
    If det = 1 Then
        Call LU_Inv(n, a(), ndx(), INV())
        If Cdet Then det = LU_Det(n, a(), ndx(), s)
    End If
    LU_FacInv = det
End Function
Function LU_FacInv2(ByVal n As Integer, _
                    ByRef a() As Double, _
           Optional ByVal Cdet As Boolean = False) As Double

    Dim i As Integer, j As Integer
    Dim s As Integer, det As Double
    Dim ndx() As Integer, INV() As Double
    
    det = LU_Fac(n, a(), ndx())
    If det = 1 Then
        Call LU_Inv(n, a(), ndx(), INV())
        If Cdet Then det = LU_Det(n, a(), ndx(), s)
        For i = 1 To n
            For j = 1 To n
                a(i, j) = INV(i, j)
            Next j
        Next i
    End If
    LU_FacInv2 = det
End Function
' ================================================
'
'            Procédures de base
'          Système à matrice pleine
'
' ================================================
Function LU_Fac(ByVal n As Integer, _
                ByRef a() As Double, _
                ByRef ndx() As Integer, _
       Optional ByRef s As Integer) As Integer

    Dim p As Integer, i As Integer, j As Integer
    Dim pp As Integer, ii As Integer
    Dim w As Double, pivot As Double
    
    If n < 2 Then Stop
    ReDim ndx(n)
    For p = 1 To n
        ndx(p) = p
    Next p
    s = 1
    For p = 1 To n - 1
        pp = ndx(p)
        For i = p + 1 To n
            ii = ndx(i)
            If Abs(a(ii, p)) > Abs(a(pp, p)) Then
                ndx(p) = ii: ndx(i) = pp: pp = ii
                s = -s
            End If
        Next i
        pivot = a(pp, p)
        If (pivot = 0) Then LU_Fac = 0: Exit Function
        For i = p + 1 To n
            ii = ndx(i)
            w = a(ii, p) / pivot
            a(ii, p) = w
            For j = p + 1 To n
                a(ii, j) = a(ii, j) - w * a(pp, j)
            Next j
        Next i
    Next p
    If a(ndx(n), n) = 0 Then LU_Fac = 0: Exit Function
    LU_Fac = 1
End Function
Function LU_Det(ByVal n As Integer, _
                ByRef LU() As Double, _
                ByRef ndx() As Integer, _
                ByVal s As Integer) As Double

    Dim p As Integer
    Dim det As Double
    
    If Abs(s) <> 1 Then Stop
    det = s
    For p = 1 To n
        det = det * LU(ndx(p), p)
    Next p
    LU_Det = det
End Function
Sub LU_Sol(ByVal n As Integer, _
           ByRef LU() As Double, _
           ByRef ndx() As Integer, _
           ByRef b() As Double, _
           ByRef x() As Double)

    Dim i As Integer, j As Integer, ii As Integer
    Dim Som As Double
    
    If n < 2 Then Stop
    ReDim x(1 To n)
    'Descente
    x(1) = b(ndx(1))
    For i = 2 To n
        ii = ndx(i)
        Som = b(ii)
        For j = 1 To i - 1
            Som = Som - LU(ii, j) * x(j)
        Next j
        x(i) = Som
    Next i
    'Remontée
    x(n) = x(n) / LU(ndx(n), n)
    For i = n - 1 To 1 Step -1
        ii = ndx(i)
        Som = x(i)
        For j = i + 1 To n
            Som = Som - LU(ii, j) * x(j)
        Next j
        x(i) = Som / LU(ii, i)
    Next i
End Sub
Sub LU_Sol2(ByVal n As Integer, _
           ByRef LU() As Double, _
           ByRef ndx() As Integer, _
           ByRef b() As Double)

    Dim i As Integer, j As Integer, ii As Integer
    Dim Som As Double, x() As Double
    
    If n < 2 Then Stop
    ReDim x(1 To n)
    'Descente
    x(1) = b(ndx(1))
    For i = 2 To n
        ii = ndx(i)
        Som = b(ii)
        For j = 1 To i - 1
            Som = Som - LU(ii, j) * x(j)
        Next j
        x(i) = Som
    Next i
    'Remontée
    b(n) = x(n) / LU(ndx(n), n)
    For i = n - 1 To 1 Step -1
        ii = ndx(i)
        Som = x(i)
        For j = i + 1 To n
            Som = Som - LU(ii, j) * b(j)
        Next j
        b(i) = Som / LU(ii, i)
    Next i
End Sub
Sub LU_Inv(ByVal n As Integer, _
           ByRef LU() As Double, _
           ByRef ndx() As Integer, _
           ByRef INV() As Double)

    Dim i As Integer, j As Integer
    Dim c As Integer, ii As Integer
    Dim Som As Double, piv As Double
    Dim x() As Double
    
    If n < 2 Then Stop
    ReDim x(1 To n), INV(1 To n, 1 To n)
    piv = LU(ndx(n), n)
    For c = 1 To n
        INV(c, c) = 1
        x(1) = INV(ndx(1), c)
        For i = 2 To n
            ii = ndx(i)
            Som = INV(ii, c)
            For j = 1 To i - 1
                Som = Som - LU(ii, j) * x(j)
            Next j
            x(i) = Som
        Next i
        x(n) = x(n) / piv
        INV(n, c) = x(n)
        For i = n - 1 To 1 Step -1
            ii = ndx(i)
            Som = x(i)
            For j = i + 1 To n
                Som = Som - LU(ii, j) * x(j)
            Next j
            x(i) = Som / LU(ii, i)
            INV(i, c) = x(i)
        Next i
    Next c
End Sub
Function LU_Residue(ByVal n As Integer, _
                    ByRef LU() As Double, _
                    ByRef b() As Double, _
                    ByRef x() As Double, _
                    ByRef r() As Double, _
           Optional ByVal Lp As Integer = 0) As Double
                          
    Dim s As Double, i As Integer, j As Integer
    
    If n < 2 Then Stop
    ReDim r(n)
    For i = n To 1 Step -1
        s = 0
        For j = i To n
            s = s + LU(i, j) * x(j)
        Next j
        r(i) = s
    Next i
    r(1) = r(1) - b(1)
    For i = 2 To n
        s = r(i)
        For j = i - 1 To 1 Step -1
            s = s + LU(i, j) * x(j)
        Next j
        r(i) = s - b(i)
    Next i
    Select Case Lp
    Case 0
        s = 0
    Case 1
        s = 0
        For i = 1 To n
            s = s + Abs(r(i))
        Next i
    Case Is < 0, Is > 10
        Stop
    Case Else
        s = 0
        For i = 1 To n
            s = s + r(i) ^ Lp
        Next i
        s = s ^ (1 / Lp)
    End Select
    LU_Residue = s
End Function
Function LU_Res(ByVal n As Integer, _
                ByRef LU() As Double, _
                ByRef b() As Double, _
                ByRef x() As Double, _
       Optional ByVal Lp As Integer = 0) As Double
                          
    Dim s As Double, i As Integer, j As Integer
    Dim r() As Double
    
    If n < 2 Then Stop
    ReDim r(n)
    For i = n To 1 Step -1
        s = 0
        For j = i To n
            s = s + LU(i, j) * x(j)
        Next j
        r(i) = s
    Next i
    r(1) = r(1) - b(1)
    For i = 2 To n
        s = r(i)
        For j = i - 1 To 1 Step -1
            s = s + LU(i, j) * x(j)
        Next j
        r(i) = s - b(i)
    Next i
    Select Case Lp
    Case 0
        s = 0
    Case 1
        s = 0
        For i = 1 To n
            s = s + Abs(r(i))
        Next i
    Case Is < 0, Is > 10
        Stop
    Case Else
        s = 0
        For i = 1 To n
            s = s + r(i) ^ Lp
        Next i
        s = s ^ (1 / Lp)
    End Select
    LU_Res = s
End Function

  


