Attribute VB_Name = "SEL"
'
'  Cours de Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module       : SEL
'  Requis       : Néant
'  Création     : Oct. 2002
'  Modification :
'  Description  : Méthodes directes de résolution
'                 d'un système linéaire Ax = b
'
Option Explicit
Function Crout_Dec(ByVal N As Integer, _
                   ByRef MAT() As Double, _
                   ByRef NDX() As Integer, _
          Optional ByRef bool As Boolean) As Integer

    Const eps = 1E-30
    Dim i As Integer, j As Integer, Pmax As Integer
    Dim Vmax As Double, v  As Double
    Dim Temp() As Double
    
    If (N <= 0) Then Stop ' Données non valides
    ReDim Temp(N)
    
    ' Calcul des éléments max de chaque ligne
    For i = 1 To N
        Vmax = 0#
        For j = 1 To N
            v = Abs(MAT(i, j))
            If (v > Vmax) Then Vmax = v
        Next j
        If Vmax = 0# Then Crout_Dec = -1: Exit Function ' matrice singulière
        Temp(i) = 1# / Vmax
    Next i
    bool = True
    ' Décomposition de Crout
    For j = 1 To N
        For i = 1 To j - 1
            Call prv_Crout_Somme(i - 1, i, j, MAT())
        Next i
        Vmax = 0#
        For i = j To N
            Call prv_Crout_Somme(j - 1, i, j, MAT())
            v = Temp(i) * Abs(MAT(i, j))
            If (v >= Vmax) Then
                Pmax = i
                Vmax = v
            End If
        Next i
        If (j <> Pmax) Then
            Call prv_Crout_Echange(N, Pmax, j, MAT())
            Temp(Pmax) = Temp(j)
            bool = Not bool
        End If
        NDX(j) = Pmax
        If (MAT(j, j) = 0#) Then MAT(j, j) = eps
        If (j <> N) Then
            v = 1# / MAT(j, j)
            For i = j + 1 To N
                MAT(i, j) = v * MAT(i, j)
            Next i
        End If
    Next j
    Crout_Dec = 1   ' OK
End Function
Sub Crout_RET(ByVal N As Integer, _
              ByRef MAT() As Double, _
              ByRef NDX() As Integer, _
              ByRef sol() As Double)

    Dim i As Integer, j As Integer
    Dim k As Integer, m As Integer
    Dim s As Double

    k = 0
    For i = 1 To N
        m = NDX(i)
        s = sol(m)
        sol(m) = sol(i)
        If (k <> 0) Then
            For j = k To i - 1
                s = s - MAT(i, j) * sol(j)
            Next j
        ElseIf (s <> 0#) Then
            k = i
        End If
        sol(i) = s
    Next i
    For i = N To 1 Step -1
        s = sol(i)
        For j = i + 1 To N
            s = s - MAT(i, j) * sol(j)
        Next j
        sol(i) = s / MAT(i, i)
    Next i
End Sub
Function Crout(ByVal N As Integer, _
               ByRef MAT() As Double, _
               ByRef b() As Double) As Integer

    Dim NDX() As Integer, err As Integer
    If N < 2 Then Stop ' Données non valides
    
    ReDim NDX(N)
    err = Crout_Dec(N, MAT(), NDX())
    If err > 0 Then Call Crout_RET(N, MAT(), NDX(), b())
    Crout = err
End Function
Function Thomas(ByVal N As Integer, _
                ByRef bas() As Double, _
                ByRef diag() As Double, _
                ByRef haut() As Double, _
                ByRef b() As Double) As Integer

    Dim gamma() As Double
    Dim den As Double
    Dim i As Integer
    
    ReDim gamma(1 To N)
    
    If N < 2 Then Stop ' Données non valides
    den = diag(1)
    If (den = 0#) Then Thomas = -1: Exit Function ' Singularité
    b(1) = b(1) / den
    For i = 2 To N
        gamma(i) = haut(i - 1) / den
        den = diag(i) - bas(i) * gamma(i)
        If (den = 0) Then Thomas = -1: Exit Function
        b(i) = (b(i) - bas(i) * b(i - 1)) / den
    Next i
    For i = N - 1 To 1 Step -1
        b(i) = b(i) - gamma(i + 1) * b(i + 1)
    Next i
    Thomas = 1 ' OK
End Function
Function Residue_Vec(ByVal N As Integer, _
                     ByRef a() As Double, _
                     ByRef b() As Double, _
                     ByRef x() As Double, _
                     ByRef r() As Double) As Double

    Dim i As Integer, j As Integer
    Dim s As Double, res As Double
    
    If N < 2 Then Stop ' Données non valides
    
    res = 0#
    For i = 1 To N
        s = 0#
        For j = 1 To N
            s = s + a(i, j) * x(j)
        Next j
        s = s - b(i)
        r(i) = s
        res = res + s * s
    Next i
    Residue_Vec = Sqr(res)
End Function
Function Residue(ByVal N As Integer, _
                 ByRef a() As Double, _
                 ByRef b() As Double, _
                 ByRef x() As Double, _
        Optional ByRef abr As Double, _
        Optional ByRef mxr As Double) As Double

    Dim i As Integer, j As Integer
    Dim s As Double, res As Double
    
    If N < 2 Then Stop ' Données non valides
    res = 0#
    abr = 0#
    mxr = 0#
    For i = 1 To N
        s = 0#
        For j = 1 To N
            s = s + a(i, j) * x(j)
        Next j
        s = Abs(s - b(i))
        res = res + s * s
        abr = abr + s
        If s > mxr Then mxr = s
    Next i
    Residue = Sqr(res)
End Function
Function Residue2(ByVal N As Integer, _
                  ByRef a() As Double, _
                  ByRef b() As Double, _
                  ByRef x() As Double) As Double

    Dim i As Integer, j As Integer
    Dim s As Double, mr As Double
    
    If N < 2 Then Stop ' Données non valides
    mr = 0#
    For i = 1 To N
        s = 0#
        For j = 1 To N
            s = s + a(i, j) * x(j)
        Next j
        s = s - b(i)
        mr = mr + s * s
    Next i
    Residue = mr
End Function
Sub GaussJordan(ByRef MAT() As Double, _
                ByVal N As Integer, _
                ByRef b() As Double, _
                ByVal m As Integer)
                                                      
    '(whp020030)
    Const SOURCE As String = "matrx::gj"
    Const ERRCODE As String = "whp020030"
    Const Nmax As Integer = 50
    Dim i    As Integer, j As Integer
    Dim k    As Integer, L As Integer
    Dim LL   As Integer
    Dim icol As Integer, irow As Integer
    Dim ndxR(1 To Nmax) As Integer
    Dim ndxC(1 To Nmax) As Integer
    Dim piv(1 To Nmax)  As Integer
    Dim big    As Double, swp As Double
    Dim pivinv As Double
    
    For j = 1 To N
        piv(j) = 0
    Next j
    For i = 1 To N
        big = 0
        For j = 1 To N
            If piv(j) <> 1 Then
                For k = 1 To N
                    If piv(k) > 1 Then Stop  ' Singularité
                    If piv(k) = 0 And Abs(MAT(j, k)) >= big Then
                            big = Abs(MAT(j, k))
                            irow = j
                            icol = k
                    End If
                Next k
            End If
        Next j
        piv(icol) = piv(icol) + 1
        If irow <> icol Then
            For L = 1 To N
                swp = MAT(irow, L):    MAT(irow, L) = MAT(icol, L):    MAT(icol, L) = swp
            Next L
            For L = 1 To m
                swp = b(irow, L):    b(irow, L) = b(icol, L):    b(icol, L) = swp
            Next L
        End If
        ndxR(i) = irow
        ndxC(i) = icol
        If MAT(icol, icol) = 0 Then Stop ' singularité
        pivinv = 1# / MAT(icol, icol)
        MAT(icol, icol) = 1#
        For L = 1 To N
            MAT(icol, L) = MAT(icol, L) * pivinv
        Next L
        For L = 1 To m
            b(icol, L) = b(icol, L) * pivinv
        Next L
        For LL = 1 To N
            If LL <> icol Then
                swp = MAT(LL, icol):   MAT(LL, icol) = 0#
                For L = 1 To N
                    MAT(LL, L) = MAT(LL, L) - MAT(icol, L) * swp
                Next L
                For L = 1 To m
                    b(LL, L) = b(LL, L) - b(icol, L) * swp
                Next L
            End If
        Next LL
    Next i
    For L = N To 1 Step -1
        If ndxR(L) <> ndxC(L) Then
            For k = 1 To N
                swp = MAT(k, ndxR(L)):   MAT(k, ndxR(L)) = MAT(k, ndxC(L)):  MAT(k, ndxC(L)) = swp
            Next k
        End If
    Next L
End Sub
' -------------------------------------------
'
'        Procédures internes de calcul
'
' -------------------------------------------
Private Sub prv_Crout_Somme(ByVal m As Integer, _
                            ByVal i As Integer, _
                            ByVal j As Integer, _
                            ByRef MAT() As Double)

    Dim k As Integer
    Dim s As Double
    
    s = MAT(i, j)
    For k = 1 To m
        s = s - MAT(i, k) * MAT(k, j)
    Next k
    MAT(i, j) = s
End Sub
Private Sub prv_Crout_Echange(ByVal N As Integer, _
                              ByVal i As Integer, _
                              ByVal j As Integer, _
                              ByRef MAT() As Double)

    Dim k    As Integer
    Dim Temp As Double
    
    For k = 1 To N
        Temp = MAT(i, k)
        MAT(i, k) = MAT(j, k)
        MAT(j, k) = Temp
    Next k
End Sub




