Attribute VB_Name = "EDP"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << EDP >>
'  Requis   :
'  Externes :
'
'  Descrip. : Equations aux Dérivées Partielles
'
Option Explicit
Function Laplace_SOR(ByVal ik1 As Integer, _
                     ByVal ik2 As Integer, _
                     ByVal jk1 As Integer, _
                     ByVal jk2 As Integer, _
                     ByRef u() As Double, _
                     ByVal eps As Double, _
                     ByVal Kmax As Integer, _
            Optional ByVal r As Double = 1#, _
            Optional ByVal w As Double = 1#, _
            Optional ByRef res1 As Double, _
            Optional ByRef res2 As Double) As Integer

    Dim i As Integer, j As Integer
    Dim N As Integer, k As Integer
    Dim bye As Boolean
    Dim r2 As Double, c As Double
    Dim anc As Double, maj As Double, dif As Double
    
    If r <= 0 Then Stop
    If w <= 0 Or w >= 2 Then Call Rho2(ik1, ik2, jk1, jk2, r, w)
    N = (ik2 - ik1 - 1) * (jk2 - jk1 - 1)
    r2 = r * r: c = 0.5 / (1 + r2): r2 = c * r2
    k = 0
    Do
        k = k + 1
        If k > Kmax Then k = -1: Exit Do
        bye = True: res1 = 0: res2 = 0
        For i = ik1 + 1 To ik2 - 1
            For j = jk1 + 1 To jk2 - 1
                anc = u(i, j)
                maj = r2 * (u(i, j + 1) + u(i, j - 1)) + c * (u(i + 1, j) + u(i - 1, j))
                dif = w * (maj - anc)
                u(i, j) = anc + dif
                res1 = res1 + Abs(dif)
                res2 = res2 + dif * dif
                If dif > eps Then bye = False
            Next j
        Next i
    Loop Until bye
    res1 = res1 / N
    res2 = Sqr(res2) / N
    Laplace_SOR = k
End Function
Function Laplace_SORD(ByVal ik1 As Integer, _
                      ByVal ik2 As Integer, _
                      ByVal jk1 As Integer, _
                      ByVal jk2 As Integer, _
                      ByRef u() As Double, _
                      ByVal eps As Double, _
                      ByVal Kmax As Integer, _
             Optional ByVal r As Double = 1#, _
             Optional ByVal w As Double = 1#, _
             Optional ByRef res1 As Double, _
             Optional ByRef res2 As Double) As Integer

    Dim i As Integer, j As Integer
    Dim N As Integer, k As Integer
    Dim bye As Boolean
    Dim r2 As Double, c As Double
    Dim anc As Double, maj As Double, dif As Double
    Dim m As Integer, s As Integer
    
    If r <= 0 Then Stop
    If w <= 0 Or w >= 2 Then Call Rho2(ik1, ik2, jk1, jk2, r, w)
    N = (ik2 - ik1 - 1) * (jk2 - jk1 - 1)
    r2 = r * r: c = 0.5 / (1 + r2): r2 = c * r2
    k = 0
    Do
        k = k + 1
        If k > Kmax Then k = -1: Exit Do
        bye = True: res1 = 0: res2 = 0
        For m = 1 To 2
            s = m
            For i = ik1 + 1 To ik2 - 1
                For j = jk1 + s To jk2 - 1 Step 2
                    anc = u(i, j)
                    maj = r2 * (u(i, j + 1) + u(i, j - 1)) + c * (u(i + 1, j) + u(i - 1, j))
                    dif = w * (maj - anc)
                    u(i, j) = anc + dif
                    res1 = res1 + Abs(dif)
                    res2 = res2 + dif * dif
                    If dif > eps Then bye = False
                Next j
                s = 3 - s
            Next i
        Next m
    Loop Until bye
    res1 = res1 / N
    res2 = Sqr(res2) / N
    Laplace_SORD = k
End Function
Function Laplace_SORT(ByVal ik1 As Integer, _
                      ByVal ik2 As Integer, _
                      ByVal jk1 As Integer, _
                      ByVal jk2 As Integer, _
                      ByRef u() As Double, _
                      ByVal eps As Double, _
                      ByVal Kmax As Integer, _
                      ByRef w As Double, _
             Optional ByVal r As Double = 1#, _
             Optional ByRef res1 As Double, _
             Optional ByRef res2 As Double) As Integer

    Dim i As Integer, j As Integer
    Dim N As Integer, k As Integer
    Dim bye As Boolean
    Dim r2 As Double, c As Double
    Dim anc As Double, maj As Double, dif As Double
    Dim m As Integer, s As Integer
    Dim fac As Double
    
    If w <= 0 Or w >= 2 Or r <= 0 Then Stop
    N = (ik2 - ik1 - 1) * (jk2 - jk1 - 1)
    r2 = r * r: c = 0.5 / (1 + r2): r2 = c * r2
    k = 0
    fac = 0.25 * Rho2(ik1, ik2, jk1, jk2, r)
    Do
        k = k + 1
        If k > Kmax Then k = -1: Exit Do
        bye = True: res1 = 0: res2 = 0
        For m = 1 To 2
            s = m
            For i = ik1 + 1 To ik2 - 1
                For j = jk1 + s To jk2 - 1 Step 2
                    anc = u(i, j)
                    maj = r2 * (u(i, j + 1) + u(i, j - 1)) + c * (u(i + 1, j) + u(i - 1, j))
                    dif = w * (maj - anc)
                    u(i, j) = anc + dif
                    res1 = res1 + Abs(dif)
                    res2 = res2 + dif * dif
                    If dif > eps Then bye = False
                Next j
                s = 3 - s
            Next i
            If w = 1 Then w = 2
            w = 1 / (1 - w * fac)
        Next m
    Loop Until bye
    res1 = res1 / N
    res2 = Sqr(res2) / N
    Laplace_SORT = k
End Function
Private Function Rho2(ByVal ik1 As Integer, _
                      ByVal ik2 As Integer, _
                      ByVal jk1 As Integer, _
                      ByVal jk2 As Integer, _
                      ByVal r As Double, _
             Optional ByRef w As Double) As Double
                     
    Dim x As Double, Y As Double, r2 As Double
    
    x = Cos(3.1415926 / (ik2 - ik1 + 1))
    Y = Cos(3.1415926 / (jk2 - jk1 + 1))
    r2 = r * r
    r2 = (x + r2 * Y) / (1 + r2)
    Rho2 = r2 * r2
    w = 2# / (1# + Sqr(1# - Rho2))
End Function

