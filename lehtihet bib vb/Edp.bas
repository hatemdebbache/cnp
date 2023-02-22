Attribute VB_Name = "EDP"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << EDP >>
'  Requis   :
'  Externes :
'
'  Descrip. : Equations aux dérivées partielles
'
Option Explicit
Function Laplace_Rec_SOR(ByVal ik1 As Integer, _
                         ByVal ik2 As Integer, _
                         ByVal jk1 As Integer, _
                         ByVal jk2 As Integer, _
                         ByRef u() As Double, _
                         ByVal eps As Double, _
                Optional ByVal w As Double = 1, _
                Optional ByVal r As Double = 1, _
                Optional ByVal kmax As Integer = 250) As Integer
                    
    Dim i As Integer, j As Integer, k As Integer
    Dim OK As Boolean
    Dim Uold As Double, Unew As Double
    Dim r2 As Double, cx As Double, cy As Double
    
    If ik1 > ik2 Or jk1 > jk2 Then Stop ' Données non valides
    If w >= 2 Or w < 1 Then Stop ' Facteur SOR incorrect
    r2 = r * r
    cx = 0.5 / (1# + r2)
    cy = r2 * cx
    k = 0
    Do
        k = k + 1
        If k > kmax Then k = 0: Exit Do
        OK = True
        For i = ik1 To ik2
            For j = jk1 To jk2
                Uold = u(i, j)
                Unew = cx * (u(i + 1, j) + u(i - 1, j)) + cy * (u(i, j + 1) + u(i, j - 1))
                Unew = Uold + w * (Unew - Uold)
                If OK Then OK = (Abs(Unew - Uold) < eps)
                u(i, j) = Unew
            Next j
        Next i
    Loop Until OK
    Laplace_Rec_SOR = k
End Function
Function Laplace_REC_Chk(ByVal ik1 As Integer, _
                         ByVal ik2 As Integer, _
                         ByVal jk1 As Integer, _
                         ByVal jk2 As Integer, _
                         ByRef u() As Double, _
                         ByVal eps As Double, _
                Optional ByVal w As Double = 1, _
                Optional ByVal r As Double = 1, _
                Optional ByVal kmax As Integer = 250) As Integer
                    
    Dim i As Integer, j As Integer, k As Integer
    Dim OK As Boolean, s As Integer
    Dim Uold As Double, Unew As Double
    Dim r2 As Double, cx As Double, cy As Double
    
    If ik1 > ik2 Or jk1 > jk2 Then Stop ' Données non valides
    If w >= 2 Or w < 1 Then Stop ' Facteur SOR incorrect
    r2 = r * r
    cx = 0.5 / (1# + r2)
    cy = r2 * cx
    k = 0
    Do
        k = k + 1
        If k > kmax Then k = 0: Exit Do
        OK = True
        s = 0
        For i = ik1 To ik2
            For j = jk1 + s To jk2 Step 2
                Uold = u(i, j)
                Unew = cx * (u(i + 1, j) + u(i - 1, j)) + cy * (u(i, j + 1) + u(i, j - 1))
                Unew = Uold + w * (Unew - Uold)
                u(i, j) = Unew
                If OK Then OK = (Abs(Unew - Uold) < eps)
            Next j
            s = 1 - s
        Next i
        s = 1
        For j = jk1 To jk2
            For i = ik1 + s To ik2 Step 2
                Uold = u(i, j)
                Unew = cx * (u(i + 1, j) + u(i - 1, j)) + cy * (u(i, j + 1) + u(i, j - 1))
                Unew = Uold + w * (Unew - Uold)
                u(i, j) = Unew
                If OK Then OK = (Abs(Unew - Uold) < eps)
            Next i
            s = 1 - s
        Next j
    Loop Until OK
    Laplace_REC_Chk = k
End Function


