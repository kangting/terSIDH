.intel_syntax noprefix

.section .rodata

.extern uint_custom_0
.extern uint_custom_1

.section .text

.global uint_custom_eq
uint_custom_eq:
    xor eax, eax
.set k, 0
.rept 25
    mov rdx, [rdi + 8*k]
    xor rdx, [rsi + 8*k]
    or rax, rdx
    .set k, k+1
.endr
    test rax, rax
    setz al
    movzx eax, al
    ret

.global uint_custom_set
uint_custom_set:
    cld
    mov rax, rsi
    stosq
    xor rax, rax
    mov ecx, 24
    rep stosq
    ret

.global uint_custom_len
uint_custom_len:
    mov ecx, 24
0:  bsr rax, [rdi + 8*rcx]
    jnz 1f
    loop 0b
    bsr rax, [rdi]
    lea rax, [rax + 1]
    cmovz rax, rcx
    ret
1:  shl ecx, 6
    lea rax, [rax + rcx + 1]
    ret

.global uint_custom_bit
uint_custom_bit:
    mov ecx, esi
    and ecx, 0x3f
    shr rsi, 6
    mov rax, [rdi + 8*rsi]
    shr rax, cl
    and rax, 1
    ret

.global uint_custom_add3
uint_custom_add3:
    mov rax, [rsi + 0]
    add rax, [rdx + 0]
    mov [rdi + 0], rax
    .set k, 1
    .rept 24
        mov rax, [rsi + 8*k]
        adc rax, [rdx + 8*k]
        mov [rdi + 8*k], rax
        .set k, k+1
    .endr
    setc al
    movzx eax, al
    ret

.global uint_custom_sub3
uint_custom_sub3:
    mov rax, [rsi + 0]
    sub rax, [rdx + 0]
    mov [rdi + 0], rax
    .set k, 1
    .rept 24
        mov rax, [rsi + 8*k]
        sbb rax, [rdx + 8*k]
        mov [rdi + 8*k], rax
        .set k, k+1
    .endr
    setc al
    movzx eax, al
    ret

.global uint_custom_random
uint_custom_random:
    test rsi, rsi
    jnz 0f
    mov esi, 200
    jmp randombytes
0:
    push r12
    push r13
    push rbx
    push rbp

    mov r12, rdi
    mov r13, rsi

    jz 1f
    mov rdi, rsi
    call uint_custom_len
    mov ebx, eax
    mov ebp, 1
    mov ecx, eax
    and ecx, 0x3f
    shl rbp, cl
    dec rbp

    mov rdi, r12
    cld
    mov ecx, 25
    xor eax, eax
    rep stosq

0:
    mov rdi, r12
    mov esi, ebx
    add esi, 7
    shr esi, 3
    call randombytes

    mov ecx, ebx
    shr ecx, 6
    cmp ecx, 25
    je 1f
    and [r12 + rcx*8], rbp
1:
    .set k, 24
    .rept 25
        mov rax, [r13 + 8*k]
        cmp [r12 + 8*k], rax
        ja 0b
        jb 2f
        .set k, k-1
    .endr
    jmp 0b
2:  pop rbp
    pop rbx
    pop r13
    pop r12
    ret

.global uint_custom_mul3_64
uint_custom_mul3_64:

    mulx r10, rax, [rsi +  0]
    mov [rdi +  0], rax

    mulx r11, rax, [rsi +  8]
    add rax, r10
    mov [rdi +  8], rax

    mulx r10, rax, [rsi + 16]
    adcx rax, r11
    mov [rdi + 16], rax

    mulx r11, rax, [rsi + 24]
    adcx rax, r10
    mov [rdi + 24], rax

    mulx r10, rax, [rsi + 32]
    adcx rax, r11
    mov [rdi + 32],rax

    mulx r11, rax, [rsi + 40]
    adcx rax, r10
    mov [rdi + 40],rax

    mulx r10, rax, [rsi + 48]
    adcx rax, r11
    mov [rdi + 48],rax

    mulx r11, rax, [rsi + 56]
    adcx rax, r10
    mov [rdi + 56],rax
    
    mulx r10, rax, [rsi +  64]
    adcx rax, r11
    mov [rdi +  64], rax

    mulx r11, rax, [rsi +  72]
    adcx rax, r10
    mov [rdi +  72], rax

    mulx r10, rax, [rsi + 80]
    adcx rax, r11
    mov [rdi + 80], rax

    mulx r11, rax, [rsi + 88]
    adcx rax, r10
    mov [rdi + 88], rax

    mulx r10, rax, [rsi + 96]
    adcx rax, r11
    mov [rdi + 96],rax

    mulx r11, rax, [rsi + 104]
    adcx rax, r10
    mov [rdi + 104],rax

    mulx r10, rax, [rsi + 112]
    adcx rax, r11
    mov [rdi + 112],rax

    mulx r11, rax, [rsi + 120]
    adcx rax, r10
    mov [rdi + 120],rax

    mulx r10, rax, [rsi +  128]
    adcx rax, r11
    mov [rdi +  128], rax

    mulx r11, rax, [rsi +  136]
    adcx rax, r10
    mov [rdi +  136], rax

    mulx r10, rax, [rsi + 144]
    adcx rax, r11
    mov [rdi + 144], rax

    mulx r11, rax, [rsi + 152]
    adcx rax, r10
    mov [rdi + 152], rax

    mulx r10, rax, [rsi + 160]
    adcx rax, r11
    mov [rdi + 160],rax

    mulx r11, rax, [rsi + 168]
    adcx rax, r10
    mov [rdi + 168],rax

    mulx r10, rax, [rsi + 176]
    adcx rax, r11
    mov [rdi + 176],rax

    mulx r11, rax, [rsi + 184]
    adcx rax, r10
    mov [rdi + 184],rax

    mulx r10, rax, [rsi +  192]
    adcx rax, r11
    mov [rdi +  192], rax

    ret


